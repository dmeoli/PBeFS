# --------------------------------------------------------------------------- #
#    conf.py                                                                  #
#                                                                             #
#    The configuration of Sphinx for the documentation of PBeFS: the          #
#    pages that docs/pages.py writes from the README, the design note and     #
#    the CHANGELOG (Markdown, by MyST), and the API reference, whose pages    #
#    docs/apiref.py writes from the XML that docs/doxygen/Doxyfile has        #
#    Doxygen extract from the headers and Breathe renders; the theme is that  #
#    of Read the Docs, with a switch to its dark version. Next to it, in the  #
#    directory doxygen, is the site of Doxygen. The workflow                  #
#    .github/workflows/docs.yml builds it with                                #
#                                                                             #
#        python3 docs/pages.py                                                #
#        sphinx-build -W -b html docs build-docs/html                         #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

import glob
import os
import re
import subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

project = 'PBeFS'
author = 'Donato Meoli'
copyright = '2026, Donato Meoli'
# the version is the last git tag, as for CMake [see cmake/DeriveVersion.cmake]
try:
    release = subprocess.run(['git', 'describe', '--tags', '--abbrev=0'],
                             cwd=ROOT, capture_output=True, text=True,
                             check=True).stdout.strip().lstrip('v')
except (OSError, subprocess.CalledProcessError):
    release = 'unreleased'
version = release

extensions = [
    'breathe',              # the Doxygen XML in Sphinx
    'myst_parser',          # the pages in Markdown
    'sphinx.ext.mathjax',   # the formulae of the comments
    'sphinx_copybutton',    # a button that copies a block of code
    'sphinx_rtd_dark_mode', # a switch between the light and the dark theme
]

# the light theme first, the dark one on the switch, which the browser
# remembers
default_dark_mode = False

# a signature longer than this goes one parameter per line, instead of
# being cut at the margin or wrapped at random
maximum_signature_line_length = 110

# the prompts of the shells are not copied with the commands
copybutton_prompt_text = r'\$ |>>> '
copybutton_prompt_is_regexp = True

source_suffix = {'.rst': 'restructuredtext', '.md': 'markdown'}
myst_enable_extensions = ['dollarmath', 'colon_fence']
myst_heading_anchors = 3
# the design note and the theory are read through their copies in
# generated/
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'branch-and-bound.md',
                    'design.md', 'theory']

# the XML that docs/apiref.py has Doxygen write, and from which it writes
# the pages of docs/api
breathe_projects = {'PBeFS': os.path.join(ROOT, 'build-docs', 'xml')}
breathe_default_project = 'PBeFS'


def mend_xml():
    """the XML of Doxygen made fit for Breathe, in place (it is a product of
    the build, which is read and never published): the attribute
    thread_local of a variable, which Doxygen writes and the schema of
    Breathe does not know, is dropped"""
    attribute = re.compile(r' thread_local="(yes|no)"')
    for path in glob.glob(os.path.join(breathe_projects['PBeFS'], '*.xml')):
        with open(path) as f:
            text = f.read()
        new = attribute.sub('', text)
        if new != text:
            with open(path, 'w') as f:
                f.write(new)


mend_xml()
primary_domain = 'cpp'
highlight_language = 'cpp'

html_theme = 'sphinx_rtd_theme'
html_title = 'Parallel Best-First Search'
# the content as wide as the signatures of the templates need, the formulae
# wider still scrolling in their box (see _static/pbefs.css); the site of
# Doxygen, which docs/doxygen/Doxyfile writes into build-docs/html/doxygen,
# is served in the directory doxygen of this one
html_static_path = ['_static']
html_css_files = ['pbefs.css']
# the sidebar stays in place while the page scrolls and scrolls on its own;
# it holds the layers and their groups, and the branch of the page down to
# its classes
html_theme_options = {'navigation_depth': 5, 'collapse_navigation': True,
                      'sticky_navigation': True}
html_context = {
    'display_github': True,
    'github_user': 'dmeoli',
    'github_repo': 'PBeFS',
    'github_version': 'main',
    'conf_py_path': '/docs/',
}

# --------------------------------------------------------------------------- #
#                              End of conf.py                                 #
# --------------------------------------------------------------------------- #
