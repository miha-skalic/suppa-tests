from scripts.calc_funcs import *
from matplotlib import pyplot

from scipy.stats.stats import pearsonr
import os

pcr_esrp1_events, pcr_ev_events = collect_pcr_psi()

savepath = '../psi/SE_fparam_esrp1/'
points = []
for my_file in os.listdir(savepath):
    if my_file.endswith('.gtf'):
        continue
    suppa1_events = load_se_events(savepath + my_file)
    points1 = make_points_pairs(pcr_esrp1_events, suppa1_events, verbose=False)
    number = (int(my_file[28:-4]))
    cor = pearsonr(*zip(*points1))[0]
    points.append((number, cor, len(points1)))




x, y, z = zip(*sorted(points))
fig, ax1 = pyplot.subplots()
ax1.plot(x, y, 'b-')
ax1.set_xlabel('Cumulative threshold')
ax1.set_ylabel('Pearson correlation', color='b')

ax2 = ax1.twinx()
ax2.plot(x, z, 'r')
ax2.set_ylabel('number of points considered', color='r')
pyplot.show()