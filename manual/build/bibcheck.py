#!/usr/bin/env python3
# --------------------------------------------------------------------------- #
#    bibcheck.py                                                              #
#                                                                             #
#    Checks that manual/chapters/Z-bibliography.md holds every work that the  #
#    repository cites, and that each of its DOIs was verified: the            #
#    references of the README (the same numbered list, with the same          #
#    anchors), the citations [[N]](...#ref-N) of the chapters and of the      #
#    pages of docs, the works that the headers cite by author and year,       #
#    e.g., "Burns et al. (JAIR 2010)", or by a doi:, and the entries          #
#    themselves (numbered 1, 2, ... without gaps, each with a DOI, an arXiv   #
#    identifier, or else a URL). A DOI or arXiv identifier is verified when   #
#    bibcheck-verified.tsv lists it with the title that Crossref (or arXiv)   #
#    gives for it and that title matches the title of the entry (a work with  #
#    none of them is listed there as none:ref-N with the reason); --verify    #
#    asks Crossref and arXiv for the ones not listed yet and adds them. Any   #
#    work missing or not verified is printed and the exit status is 1.        #
#                                                                             #
#        python3 manual/build/bibcheck.py [--verify]                          #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Università di Pisa                              #
# --------------------------------------------------------------------------- #

import difflib
import glob
import json
import os
import re
import sys
import unicodedata
import urllib.parse
import urllib.request

BUILD = os.path.dirname(os.path.abspath(__file__))
MANUAL = os.path.dirname(BUILD)
ROOT = os.path.dirname(MANUAL)
BIB = os.path.join(MANUAL, 'chapters', 'Z-bibliography.md')
VERIFIED = os.path.join(BUILD, 'bibcheck-verified.tsv')
AGENT = 'bibcheck/1.0 (bibliography check of a manual)'

ENTRY = re.compile(r'^<a id="ref-(\d+)"></a>\[(\d+)\] (.*?)$', re.M)
DOI = re.compile(r'(?:doi\.org/|doi:\s*)(10\.\d{4,9}/[^\s\]>"]+)', re.I)
ARXIV = re.compile(r'arxiv\.org/abs/(\d{4}\.\d{4,5})|arXiv:(\d{4}\.\d{4,5})')
URL = re.compile(r'https?://[^\s)\]>]+')
# a citation by author and year in the comments of a header, e.g., "Korf
# (1985)", "Hart, Nilsson and Raphael (1968)", "Burns et al. (JAIR 2010)"
NAME = r"[A-Z][A-Za-zÀ-ÿ'\-]+"
AUTHORYEAR = re.compile(
    r'(%s)(?:, %s)*(?: et al\.| and %s)? \((?:[^()]*?[ ,])?((?:19|20)\d\d)\)'
    % (NAME, NAME, NAME))
# a citation inside parentheses, e.g., "(Rihani, Sanders and Dementiev,
# 2015)"
PARENYEAR = re.compile(r'\((%s)(?:, %s)*(?: and %s)?, ((?:19|20)\d\d)\)'
                       % (NAME, NAME, NAME))


def doi(d):
    """a DOI matched in a text without the punctuation after it, i.e., the
    closing parentheses that it does not open and a final . , or ;"""
    d = d.rstrip('.,;')
    while d.endswith(')') and d.count(')') > d.count('('):
        d = d[:-1].rstrip('.,;')
    return d.lower()


def plain(s):
    """a string without accents, markup and punctuation, in lower case"""
    s = unicodedata.normalize('NFKD', s)
    s = ''.join(c for c in s if not unicodedata.combining(c))
    return re.sub(r'[^a-z0-9]+', ' ', s.lower()).strip()


def entries():
    """{ number : ( first author's surname , year , title , text ) }"""
    text = open(BIB, encoding='utf-8').read()
    out = {}
    for m in ENTRY.finditer(text):
        n, body = int(m.group(1)), m.group(3)
        if int(m.group(2)) != n:
            sys.exit('bibcheck: the anchor ref-%s is on the entry [%s]'
                     % (m.group(1), m.group(2)))
        if n in out:
            sys.exit('bibcheck: the entry [%d] is there twice' % n)
        authors = body.split(', [')[0].split(', *')[0]
        first = authors.split(' and ')[0].split(', ')[0].split()[-1]
        years = re.findall(r'\b((?:19|20)\d\d)\b', body)
        title = re.search(r'\*(.+?)\*', body)
        out[n] = (plain(first), years[-1] if years else '',
                  title.group(1) if title else '', body)
    return out


def verified():
    """{ identifier : title } of bibcheck-verified.tsv"""
    out = {}
    if os.path.exists(VERIFIED):
        for line in open(VERIFIED, encoding='utf-8'):
            if line.strip() and not line.startswith('#'):
                key, title = line.rstrip('\n').split('\t', 1)
                out[key.lower()] = title
    return out


def ask(key):
    """the title that Crossref (a DOI) or arXiv (arxiv:<id>) gives"""
    if key.startswith('arxiv:'):
        url = ('http://export.arxiv.org/api/query?id_list='
               + key.split(':', 1)[1])
        req = urllib.request.Request(url, headers={'User-Agent': AGENT})
        xml = urllib.request.urlopen(req, timeout=30).read().decode()
        titles = re.findall(r'<title>(.*?)</title>', xml, re.S)
        return ' '.join(titles[-1].split()) if len(titles) > 1 else None
    url = 'https://api.crossref.org/works/' + urllib.parse.quote(key)
    req = urllib.request.Request(url, headers={'User-Agent': AGENT})
    msg = json.load(urllib.request.urlopen(req, timeout=30))['message']
    return ' '.join(msg['title'][0].split()) if msg.get('title') else None


def similar(a, b):
    """whether two titles are the same work (case, accents, punctuation
    and a subtitle left out do not count)"""
    a, b = plain(a), plain(b)
    return (a == b or a.startswith(b) or b.startswith(a)
            or difflib.SequenceMatcher(None, a, b).ratio() > 0.85)


def main():
    verify = '--verify' in sys.argv[1:]
    bib = entries()
    bad = []
    nums = sorted(bib)
    if nums != list(range(1, len(nums) + 1)):
        bad.append('the entries are not numbered 1 to %d without gaps'
                   % len(nums))

    # the README holds the same list
    readme = open(os.path.join(ROOT, 'README.md'), encoding='utf-8').read()
    for m in ENTRY.finditer(readme):
        n = int(m.group(1))
        if n not in bib:
            bad.append('README [%d] is not in the bibliography' % n)
        elif m.group(3).strip() != bib[n][3].strip():
            bad.append('README [%d] differs from the bibliography' % n)
    if set(int(m.group(1)) for m in ENTRY.finditer(readme)) != set(bib):
        bad.append('the README and the bibliography list other works')

    # the citations by number of the manual and of the pages of docs
    sources = (glob.glob(os.path.join(MANUAL, 'chapters', '*.md'))
               + [os.path.join(MANUAL, 'index.md')]
               + glob.glob(os.path.join(ROOT, 'docs', '**', '*.md'),
                           recursive=True))
    for path in sources:
        text = open(path, encoding='utf-8').read()
        for m in re.finditer(r'#ref-([0-9A-Za-z_-]+)\)', text):
            key = m.group(1)
            if not key.isdigit() or int(key) not in bib:
                bad.append('%s cites ref-%s, which is not in the '
                           'bibliography' % (os.path.relpath(path, ROOT),
                                             key))

    # the works that the headers cite by author and year, or by doi
    index = {(a, y) for a, y, _, _ in bib.values()}
    dois = {doi(d) for e in bib.values() for d in DOI.findall(e[3])}
    for path in sorted(glob.glob(os.path.join(ROOT, 'include', '**', '*.h'),
                                 recursive=True)):
        text = open(path, encoding='utf-8').read()
        text = re.sub(r'\s*\n\s*\*?\s*', ' ', text)
        cites = [m for m in AUTHORYEAR.finditer(text)
                 if not m.group(1).isupper()]
        cites += list(PARENYEAR.finditer(text))
        for m in cites:
            if (plain(m.group(1)), m.group(2)) not in index:
                bad.append('%s cites %s, which is not in the bibliography'
                           % (os.path.relpath(path, ROOT), m.group(0)))
        for d in DOI.findall(text):
            if doi(d) not in dois:
                bad.append('%s cites doi:%s, which is not in the '
                           'bibliography' % (os.path.relpath(path, ROOT), d))

    # each entry with a verified DOI or arXiv identifier, or else a URL
    known = verified()
    added = []
    nodoi = []
    for n in nums:
        _, _, title, body = bib[n]
        keys = [doi(d) for d in DOI.findall(body)]
        keys += ['arxiv:' + (a or b) for a, b in ARXIV.findall(body)]
        if not keys:
            nodoi.append(n)
            if not URL.search(body) and 'none:ref-%d' % n not in known:
                bad.append('[%d] has neither a DOI, nor an arXiv '
                           'identifier, nor a URL, and bibcheck-verified.tsv '
                           'does not say why (none:ref-%d)' % (n, n))
            continue
        for key in keys:
            if key not in known and verify:
                try:
                    t = ask(key)
                except Exception as e:  # noqa: BLE001 (network errors)
                    bad.append('[%d] %s could not be verified: %s'
                               % (n, key, e))
                    continue
                if t is None:
                    bad.append('[%d] %s is unknown to its registry'
                               % (n, key))
                    continue
                known[key] = t
                added.append((key, t))
            if key not in known:
                bad.append('[%d] %s is not verified (run with --verify)'
                           % (n, key))
            elif title and not similar(title, known[key]) and not any(
                    similar(known[key], plain(x)) for x in body.split('*')):
                bad.append('[%d] %s is "%s" for its registry, not "%s"'
                           % (n, key, known[key], title))
    if added:
        new = not os.path.exists(VERIFIED)
        with open(VERIFIED, 'a', encoding='utf-8') as f:
            if new:
                f.write('# identifier\ttitle given by Crossref or arXiv\n')
            for key, t in added:
                f.write('%s\t%s\n' % (key, t))

    for b in bad:
        print('bibcheck: ' + b)
    print('bibcheck: %d entries, %d without DOI or arXiv identifier (%s), '
          '%d problems' % (len(nums), len(nodoi),
                           ', '.join(map(str, nodoi)) or 'none', len(bad)))
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()

# --------------------------------------------------------------------------- #
# ----------------------------- END bibcheck.py ----------------------------- #
# --------------------------------------------------------------------------- #
