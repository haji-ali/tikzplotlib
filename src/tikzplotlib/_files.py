from pathlib import Path
import warnings
import numpy as np

class DataFile(object):
    def __init__(self, tablename, filename, transpose=False, sep=',',
                 fillvalue="", load=False, column_fmt=None):
        self.columns = dict()
        self.tablename = tablename
        self.filename = filename
        self.transpose = transpose
        self.sep = sep
        self.fillvalue = fillvalue
        self.column_fmt = None
        if load:
            self.load()

    def append(self, col_type, column, rel_tol=1e-09, allow_partial=True):
        cmp_eq = lambda a, b: np.abs(a-b) <= rel_tol * np.maximum(np.abs(a), np.abs(b))

        ac = np.array(column)
        for j, v in self.columns.items():
            m = np.minimum(len(v), len(ac))
            if (allow_partial or len(v) == len(column)) and \
               np.all(cmp_eq(v[:m], ac[:m])):
                if len(ac) > len(v):
                    self.columns[j] = ac
                return j

        key = len(self.columns)
        if self.column_fmt is not None:
            ind = 0
            while True:
                newkey = self.column_fmt.format(ind, type=col_type, ind=ind)
                if key == newkey:
                    raise ValueError("DataFile.column_fmt must contain '{ind}' ")
                key = newkey
                if key not in self.columns:
                    break
                ind += 1
        assert(key not in self.columns)
        # New column
        self.columns[key] = ac
        return key

    def load(self):
        ### TODO: Does not work with fillvalue other than nan
        ### TODO: Does not work non-labelled data
        if not os.path.isfile(self.filename):
            return
        raise Exception("Not updated")
        with open(self.filename, 'r') as f:
            if self.transpose:
                for line in f:
                    parts = line.strip().split(self.sep)
                    while parts[-1] == self.fillvalue:
                        parts.pop()
                    self.columns[parts[0]] = np.array(parts[1:]).astype(np.float)
            else:
                keys = f.readline().strip().split(self.sep)
                values = []
                for line in f:
                    values = line.strip().split(self.sep)
                    values.append(np.array(values).astype(np.float))
                values = np.vstack(values)
                for i, k in enumerate(keys):
                    l = len(values[:, k])
                    while values[l-1, k] == self.fillvalue:
                        l -= 1
                    self.columns[k] = values[:l, k]

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.write()

    def write(self):
        if len(self.columns) == 0:
            warnings.warn("Datafile is empty")
            return

        if not self.transpose:
            # TODO: Implement this
            with open(self.filename, 'w') as f:
                for row in zip_longest(*self.columns.values(), fillvalue=self.fillvalue):
                    f.write(self.sep.join(
                        ["%.15g" % v if isinstance(v, numbers.Number)
                         else str(v) for v in row]) + "\n")
        else:
            max_count = np.max([len(v_row) for v_row in self.columns.values()])
            with open(self.filename, 'w') as f:
                for v_label, v_row in self.columns.items():
                    data = []
                    if self.column_fmt is not None:
                        data.append(v_label)
                    data.extend(["%.15g" % v for v in v_row])

                    # Pad with fill-values... Pgfplot needs it for some reason
                    if len(v_row) < max_count:
                        data.extend([str(self.fillvalue)] * (max_count - len(v_row)))
                    f.write(self.sep.join(data))
                    f.write("\n")


def _gen_filepath(data, nb_key, ext):
    rel_filepath = Path(f"{data['base name']}-{data[nb_key]:03d}{ext}")

    if data["rel data path"]:
        rel_filepath = data["rel data path"] / rel_filepath

    return data["output dir"] / rel_filepath, rel_filepath


def new_filepath(data, file_kind, ext):
    """Returns an available filepath.

    :param file_kind: Name under which numbering is recorded, such as 'img' or
                      'table'.
    :type file_kind: str

    :param ext: Filename extension.
    :type ext: str

    :returns: (filepath, rel_filepath) where filepath is a path in the
              filesystem and rel_filepath is the path to be used in the tex
              code.
    """

    nb_key = file_kind + "number"
    if nb_key not in data.keys():
        data[nb_key] = -1

    data[nb_key] += 1
    filepath, rel_filepath = _gen_filepath(data, nb_key, ext)
    if not data["override externals"]:
        # Make sure not to overwrite anything.
        while filepath.is_file():
            data[nb_key] += 1
            filepath, rel_filepath = _gen_filepath(data, nb_key, ext)

    return filepath, rel_filepath
