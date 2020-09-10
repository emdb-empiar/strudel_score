from math import fabs
import time


def find_y_label_coordinates(lst, y_range, widget_height, label_size, spacing=0.9):
    lst = [i for i in lst if y_range[0] < i < y_range[1]]
    start = time.time()
    if len(lst) > 0:
        while label_size > 2:
            text_h = (y_range[1] - y_range[0]) / widget_height * label_size
            line_height = (spacing + 1) * text_h
            result = spread_points(lst, line_height, y_range)

            label_size -= 0.5
            if result is not None:
                print('time', time.time()-start)
                return result, label_size
        else:
            return lst, label_size


def spread_points(lst, line_height, y_range):
    H = line_height
    total_height = (len(lst)-1) * H
    if total_height > y_range[1] - y_range[0]:
        return None
    err = line_height / 10

    T = y_range[1] - H
    B = y_range[0] + H
    overlaps = True
    z = 0
    while overlaps:
        overlaps = False
        for i, val in enumerate(lst[:-1]):
            z += 1
            d = val - lst[i+1]
            if d + err < H:
                overlaps = True
                shift = (H - d) / 2

                if i == 0 and T - val < shift:
                    lst[i] = T
                elif i != 0 and lst[i-1] - lst[i] < shift:
                    lst[i] = (lst[i - 1] + lst[i]) / 2
                else:
                    lst[i] = val + shift

                if i == len(lst)-2 and val - B < shift:
                    lst[i+1] = B
                elif i != len(lst)-2 and lst[i+1] - lst[i+2] < shift:
                    lst[i + 1] = (lst[i + 2] + lst[i + 1]) / 2
                else:
                    lst[i + 1] = lst[i + 1] - shift
    return lst


def cluster_spreed_points(lst, min_d):
    """
    Takes a list of numbers and makes the distance between any adjacent
    number greater than min_d. If the output values form clusters
    :param lst:
    :param min_d:
    :return:
    """
    # Create clusters of points
    clusters = []
    cluster = [lst[0]]
    for i, val in enumerate(lst[:-1]):
        if fabs(val - lst[i + 1]) < min_d:
            cluster.append(lst[i + 1])
        else:
            clusters.append(cluster)
            cluster = [lst[i + 1]]
    clusters.append(cluster)
    print(clusters)
    # Merge clusters which overlap
    if len(clusters) > 2:
        tmp = clusters
        while True:
            merged = []
            for i, c1 in enumerate(tmp[:-1]):
                c2 = tmp[i + 1]
                d_c1 = min(c1) + (max(c1) - min(c1)) / 2
                d_c2 = min(c2) + (max(c2) - min(c2)) / 2
                d = fabs(d_c1 - d_c2)

                if d < (len(c1) / 2 + len(c2) / 2) * min_d:
                    merged.append(c1 + c2)
                    merged = merged + tmp[i + 2:]
                    break
                else:
                    merged.append(c1)
            else:
                merged.append(tmp[-1])
            if merged == tmp:
                break
            else:
                tmp = merged
    else:
        merged = clusters

    # Calculate label positions
    out = []
    for cluster in merged:
        ln = len(cluster)
        if ln > 1:
            mdl_v = min(cluster) + (max(cluster) - min(cluster)) / 2
            for i, v in enumerate(cluster):
                if v > mdl_v:
                    val = mdl_v + (ln / 2 - i - 0.5) * min_d
                else:
                    val = mdl_v - (i - ln / 2 + 0.5) * min_d
                out.append(val)
        else:
            out.append(cluster[0])
    return out
