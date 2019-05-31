import pstats
from pstats import SortKey
p = pstats.Stats('../log/190529process_isoloci1.profile')
p.sort_stats(SortKey.CUMULATIVE).print_stats(25)
