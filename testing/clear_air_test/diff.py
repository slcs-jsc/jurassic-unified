import sys
import glob
import re

def disassemble(lines):
    comments, count = -1, 0
    ret = []
    for l in lines:
        if len(l) > 1:
            if l[0] == '#':
                count += 1
            else:
                if comments == -1:
                    comments = count
                assert comments == count
                assert comments % 2 == 0
                row = list(map(float, l.strip().split()))
                # assert len(row) == comments
                separated = [row[:10], row[10:]]
                ret.append(separated)
    return ret

def calc_diff(mat, a, b):
    assert len(a) >= len(b)
    assert len(a) % 2 == 0 and len(b) % 2 == 0
    half_a = int(len(a) / 2)
    half_b = int(len(b) / 2)
    max_abs_diff = 0
    max_rel_diff = 0
    pairs = []
    for i in range(0, len(b)): #check both tau and rad
        id_a = i
        if i >= half_b:
            id_a = half_a + i - half_b
        abs_diff = abs(a[id_a] - b[i])
        pairs.append((a[id_a], b[i]))
        if abs(a[id_a]) < 1e-15:
            rel_diff = 1e9
        else:
            rel_diff = abs((b[i] - a[id_a]) / a[id_a])
        if abs_diff > max_abs_diff:
            max_abs_diff = abs_diff
        if rel_diff > max_rel_diff:
            max_rel_diff = rel_diff
    mat.append(pairs)
    return [max_abs_diff, max_rel_diff]

#a and b are file paths of outputs we want to compare
def check(a, b):
    print("comparing {} and {}".format(a, b))
    mat = []
    assert a != b, "arguments must differ!"
    lines_a, lines_b = None, None
    try:
        with open(a) as f:
            lines_a = f.readlines()
    except IOError:
        raise Exception("First file not accessible")

    try:
        with open(b) as f:
            lines_b = f.readlines()
    except IOError:
        raise Exception("Second file not accessible")

    d_a = disassemble(lines_a)
    d_b = disassemble(lines_b)

    max_abs_diff, max_rel_diff = 0, 0
    max_input_diff = 0
    assert len(d_a) >= len(d_b)
    for i in range(len(d_b)):
        assert len(d_a[i][0]) == len(d_b[i][0])
        for j in range(len(d_a[i][0])):
            t = abs(d_a[i][0][j] - d_b[i][0][j])
            if t > max_input_diff:
                max_input_diff = t
        #assert d_a[i][0] == d_b[i][0]
        abs_diff, rel_diff = calc_diff(mat, d_a[i][1], d_b[i][1])
        if abs_diff > max_abs_diff:
            max_abs_diff = abs_diff
        if rel_diff > max_rel_diff:
            max_rel_diff = rel_diff
    print("number of observations:", len(d_b))
    print("number of channels:", int(len(d_b[0][1]) / 2))
    print("max input diff: {}".format(max_input_diff));
    print("max abs diff: {}, max rel dif: {}".format(max_abs_diff, max_rel_diff))
    if len(mat) > 0:
      print("mat size: {} x {}".format(len(mat), len(mat[0])));
      for i in range(min(len(mat), 5)):
        for j in range(min(len(mat[0]), 5)):
          print(mat[i][j], end=" ")
        print()

if __name__ == "__main__":
    assert len(sys.argv) == 2
    print('check("rad-785-812-after-bugfix.org, {}")'.format(sys.argv[1]))
    check("rad-785-812-after-bugfix.org", sys.argv[1])
