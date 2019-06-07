import pstats
from pstats import SortKey
p = pstats.Stats('../log/190605tabu.profile')
p.sort_stats(SortKey.CUMULATIVE).print_stats(25)
