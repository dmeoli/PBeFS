###############################################################################
#                                                                             #
#   ratios.py: the ratios of times as the tables of the README write them     #
#                                                                             #
#   A ratio r of two times, that of a run over that of the reference, is      #
#   written as how many times the reference is faster ("2.0× faster") or      #
#   slower ("1.6× slower"), and as "tie" when it is within the A/A noise      #
#   NOISE of 1; a range of ratios as "1.2× to 3.4× faster", "1.1× slower to   #
#   14× faster" or "tie to 2.1× faster"; a speedup as "3.8×". A number has    #
#   one decimal below 9.95 and none from there on, so that 1.95 is 2.0 and    #
#   9.96 is 10. With html, the sign × is the entity &times;.                  #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
#                                                                             #
###############################################################################

NOISE = 0.05


def times(x, html=False):
    """the number x of times, as "2.0×" or "14×\""""
    x += 1e-9
    return ('%.0f' % x if x >= 9.95 else '%.1f' % x) + \
        ('&times;' if html else '×')


def verdict(r, html=False):
    """the ratio r of the time of a run over that of the reference"""
    if abs(r - 1) <= NOISE:
        return 'tie'
    if r > 1:
        return times(r, html) + ' faster'
    return times(1 / r, html) + ' slower'


def span(lo, hi, html=False):
    """the range [ lo , hi ] of the ratios of the time of a run over that of
    the reference"""
    a, b = verdict(lo, html), verdict(hi, html)
    if a == b:
        return a
    if lo > 1 + NOISE:                       # faster at both ends
        return '%s to %s faster' % (times(lo, html), times(hi, html))
    if hi < 1 - NOISE:                       # slower at both ends
        return '%s to %s slower' % (times(1 / hi, html), times(1 / lo, html))
    return '%s to %s' % (a, b)

###############################################################################
############################## End of ratios.py ###############################
###############################################################################
