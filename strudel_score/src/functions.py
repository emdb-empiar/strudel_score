from math import fabs
import time
import requests
import os
import re


CHUNK_SIZE = 512 * 1024

def find_y_label_coordinates(lst, y_range, widget_height, label_size, spacing=0.9):
    lst = [i for i in lst if y_range[0] < i < y_range[1]]
    if len(lst) > 0:
        while label_size > 2:
            text_h = (y_range[1] - y_range[0]) / widget_height * label_size
            line_height = (spacing + 1) * text_h
            result = spread_points(lst, line_height, y_range)

            label_size -= 0.5
            if result is not None:
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


def g_drive_download(url, out_path):
    orig_url = url
    sess = requests.session()
    headers = {
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/39.0.2171.95 Safari/537.36"}
    while True:
        res = sess.get(url, headers=headers, stream=True)
        if "Content-Disposition" in res.headers:
            # This is the file
            break
        url = get_url_from_gdrive_confirmation(res.text)

        if url is None:
            print("Permission denied:", orig_url)
            print(
                "Maybe you need to change permission over "
                "'Anyone with the link'?")
            return

    with open(out_path, "wb") as f:
        # total = res.headers.get("Content-Length")
        size = 0
        for chunk in res.iter_content(chunk_size=CHUNK_SIZE):
            if total is not None:
                total = int(total)
            f.write(chunk)
            print((size + len(chunk) / 1000000))
            size += len(chunk)


def get_url_from_gdrive_confirmation(contents):
    url = ""
    for line in contents.splitlines():
        m = re.search(r'href="(\/uc\?export=download[^"]+)', line)
        if m:
            url = "https://docs.google.com" + m.groups()[0]
            url = url.replace("&amp;", "&")
            return url
        m = re.search("confirm=([^;&]+)", line)
        if m:
            confirm = m.groups()[0]
            url = re.sub(
                r"confirm=([^;&]+)", r"confirm={}".format(confirm), url
            )
            return url
        m = re.search('"downloadUrl":"([^"]+)', line)
        if m:
            url = m.groups()[0]
            url = url.replace("\\u003d", "=")
            url = url.replace("\\u0026", "&")
            return url
        m = re.search('<p class="uc-error-subcaption">(.*)</p>', line)
        if m:
            error = m.groups()[0]
            raise RuntimeError(error)


def unpack_tar(file, output):
    import tarfile
    my_tar = tarfile.open(os.path.abspath(file))
    my_tar.extractall(output)
    my_tar.close()


def read_gfile(url):
    orig_url = url
    sess = requests.session()
    headers = {
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/39.0.2171.95 Safari/537.36"}
    while True:
        res = sess.get(url, headers=headers, stream=True)
        if "Content-Disposition" in res.headers:
            # This is the file
            break
        url = get_url_from_gdrive_confirmation(res.text)

        if url is None:
            print("Permission denied:", orig_url)
            print(
                "Maybe you need to change permission over "
                "'Anyone with the link'?")
            return

    chunks = [chunk for chunk in res.iter_content(chunk_size=CHUNK_SIZE)]
    return chunks[0].decode("utf-8")


def read_config_value(path, key):
    try:
        with open(path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                words = line.partition(' = ')
                if key == words[0]:
                    return words[-1]
    except FileNotFoundError:
        return None



def set_library(lib_path):

    lib_path = os.path.abspath(os.path.expanduser(os.path.expandvars(lib_path)))
    if not os.path.exists(lib_path):
        raise Exception(f'Path {lib_path} does not exist!')

    files = os.listdir(lib_path)
    files = [f for f in files if f.startswith('motifs_')]
    if len(files) == 0:
        raise Exception("Error", f"Could not find strudel libraries in the directory {lib_path}\n"
                         f"Please select a directory which contain folders named motifs_***")

    record = f'lib_path={lib_path}\n'
    lines = [record]
    basedir = os.path.dirname(os.path.abspath(__file__))
    config_path = os.path.join(basedir, 'config.txt')
    try:
        with open(config_path, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if 'lib_path' in line:
                    lines[i] = record
    except FileNotFoundError:
        pass
    with open(config_path, 'w') as f:
        for line in lines:
            f.write(line)
    print('Info', f'Strudel library set to: {lib_path}')


