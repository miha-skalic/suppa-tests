from scipy.stats.stats import pearsonr
from matplotlib import pyplot


def load_sepcr_events(file_path):
    """
    Loads in events from file ('the new format')
    """
    event_dict = {}

    with open(file_path) as handle:
        handle.__next__()
        for line in handle:
            line = line.strip().split('\t')
            exp_val = line[-1]
            line = line[0].split(';')[1]
            locs = line.split(':')
            position = (locs[1], int(locs[2]), int(locs[3]))
            event_dict[position] = event_dict.get(position, []) + [float(exp_val)]
    return event_dict

def load_se_events(file_path):
    """
    Loads in events from file
    """
    event_dict = {}

    with open(file_path) as handle:
        handle.__next__()
        for line in handle:
            line = line.strip().split('\t')
            exp_val = line[-1]
            line = line[0].split(';')[1]
            position = tuple([line.split(':')[1]] + [int(x) for x in line.split('-')[1].split(':')])
            if exp_val.startswith('NA'):
                continue
            event_dict[position] = event_dict.get(position, []) + [float(exp_val)]
    return event_dict


def collect_pcr_psi():
    pcr_ev = {}
    pcr_esrp1 = {}
    with open('../psi/pcr_psi.txt') as handle:
        for line in handle:
            line = line.strip()
            line = line.split('\t')
            flanks = line[2].split(' ')[0], int(line[2].split(' ')[2]) + 1, int(line[2].split(' ')[4])
            pcr_esrp1[flanks] = float(line[-2])
            pcr_ev[flanks] = float(line[-1])
    return pcr_esrp1, pcr_ev


def make_points_pairs(pcr_data, calculation_data, verbose=True, table_file=''):
    """
    Returns pairs of points
    """
    n_multi = 0
    point_pairs = []
    
    if table_file:
        table_file = open(table_file, 'w')
        print('Event_coordinates\tNGS_PSI\tPCR_PSI', file=table_file)
    
    for event in pcr_data:
        if event not in calculation_data:
            continue
        elif len(calculation_data[event]) == 1:
            # skip if it not expressed
            if calculation_data[event][0] < 0:
                continue
            point_pairs.append((calculation_data[event][0], pcr_data[event]))
            if table_file:
                print('{}\t{}\t{}'.format('{}:{}-{}'.format(event[0], event[1], event[2]), calculation_data[event][0], pcr_data[event]), file=table_file)
        else:
            n_multi += 1
    if verbose:
        print(('\tPoints covered: {} ({:.4f} %; removed {}'
               ' point(s) duo to "duplicated" event)').format(len(point_pairs), len(point_pairs) * 100 / len(pcr_data),
                                                              n_multi))
    return point_pairs


def plot_corelation(splot, group1, group2=None, sam1_lab='Sample1', sam2_lab='Sample2',
                    x_lab='calculated PSI', title='x'):
    """
    Produces realtions plot
    """
    splot.set_title(title)
    splot.plot([0, 1], [0, 1], 'black', alpha=0.3)
    x, y = zip(*group1)
    # plot points
    splot.plot(x, y, 'ob', label='{} - P.corr.: {:.4f}'.format(sam1_lab, pearsonr(x, y)[0]))

    if group2:
        x, y = zip(*group2)
        # plot points
        splot.plot(x, y, 'xr', label='{} - P.corr.: {:.4f}'.format(sam2_lab, pearsonr(x, y)[0]),
                   mew=1, ms=12)

    splot.set_xlabel(x_lab)
    splot.set_ylabel('PCR PSI')
    splot.set_ylim(-0.1, 1.1)
    splot.set_xlim([-0.1, 1.1])
    splot.legend(loc='best')


if __name__ == '__main__':

    # makes setup
    my_figure = pyplot.figure()
    axes = my_figure.add_subplot(121)
    axes2 = my_figure.add_subplot(122)

    pcr_esrp1_events, pcr_ev_events = collect_pcr_psi()

    print('ESRP1 SAMPLE:')
    print('SUPA v1.2: PSI from sailfish numbers:')
    suppa1_events = load_se_events('../psi/refseq_sailfish_suppa1_esrp1.txt.psi')
    points1 = make_points_pairs(pcr_esrp1_events, suppa1_events)
    print('SUPA SEpcr: PSI from sailfish numbers:')
    suppa2_events = load_sepcr_events('../psi/refseq_sailfish_SEpcr_esrp1.txt.psi')
    points2 = make_points_pairs(pcr_esrp1_events, suppa2_events)
    print()

    print('EV SAMPLE:')
    print('SUPA v1.2: PSI from sailfish numbers:')
    suppa1_events = load_se_events('../psi/refseq_sailfish_suppa1_ev.txt.psi')
    points3 = make_points_pairs(pcr_ev_events, suppa1_events)
    print('SUPA SEpcr: PSI from sailfish numbers:')
    suppa2_events = load_sepcr_events('../psi/refseq_sailfish_SEpcr_ev.txt.psi')
    points4 = make_points_pairs(pcr_ev_events, suppa2_events)

    plot_corelation(axes, points1, points2, title='ESRP1 data on Sailfish',
                    sam1_lab='suppa_1.1', sam2_lab='suppa_1.2')
    plot_corelation(axes2, points3, points4, title='EV data on Sailfish',
                    sam1_lab='suppa_1.1', sam2_lab='suppa_1.2')

    pyplot.show()