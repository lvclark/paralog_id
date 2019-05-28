import pstats
from pstats import SortKey
p = pstats.Stats('../log/process_isoloci_profile.txt')
p.sort_stats(SortKey.CUMULATIVE).print_stats(10)
