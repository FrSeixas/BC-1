import argparse
import numpy as np
from Bio.SubsMat.MatrixInfo import blosum50


class SmithWaterman:
    grid = None
    high = 0
    high_pos = None
    traces = None

    sequence_a = None
    sequence_b = None
    gap_pen = None

    def __init__(self, gap_pen, sequence_a, sequence_b) -> None:
        self.sequence_a = sequence_a
        self.sequence_b = sequence_b
        self.gap_pen = gap_pen

        self.grid = np.empty(
            (len(sequence_a) + 1, len(sequence_b) + 1), dtype=GridCell)

        for i in range(len(sequence_a) + 1):
            self.grid[i][0] = GridCell(0)

        for i in range(len(sequence_b) + 1):
            self.grid[0][i] = GridCell(0)

    def search_blosum50(self, char1, char2):
        pair = (char1, char2)

        if pair not in blosum50:
            return blosum50[(tuple(reversed(pair)))]
        return blosum50[pair]

    def solve(self):
        for j in range(1, len(self.sequence_b) + 1):
            for i in range(1, len(self.sequence_a) + 1):
                surrounding = []

                surrounding.append(((i-1, j-1), self.grid[i-1][j-1].value + self.search_blosum50(
                    self.sequence_a[i-1], self.sequence_b[j-1])))

                surrounding.append(
                    ((i-1, j), self.grid[i-1][j].value - self.gap_pen))

                surrounding.append(
                    ((i, j-1), self.grid[i][j-1].value - self.gap_pen))

                value = 0

                for cell in surrounding:
                    value = max(value, cell[1])

                arrows = [i for i, j in surrounding if j == value]

                self.grid[i][j] = GridCell(
                    value, arrows) if value else GridCell()

                if value > self.high:
                    self.high_pos = [(i, j)]
                    self.high = value

                elif value == self.high and value != 0:
                    self.high_pos.append((i, j))

    def traceback(self):
        if self.traces:
            return self.traces

        results = list()

        for pos in self.high_pos:
            results += self.traceback_from(pos)

        self.traces = results

        return self.traces

    def traceback_from(self, pos):
        traces = []

        if not pos:
            return []

        cell = self.grid[pos[0]][pos[1]]

        if cell.arrows:
            for previous in cell.arrows:
                recurs_traces = self.traceback_from(previous)
                if recurs_traces:
                    traces += [[pos] + i for i in recurs_traces]

            return traces

        else:
            return [[pos]]


class GridCell:
    arrows = None
    value = None

    def __init__(self, value=0, arrows=None) -> None:
        self.arrows = arrows
        self.value = value

    def __repr__(self):
        return f"{self.arrows}"


def read_file_content(file):
    if not file:
        raise FileNotFoundError
    with open(file, "r") as f:
        return f.readline()


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate best local alignment between two Sequences using the Smith-Waterman Algorithm.")

    parser.add_argument("-p", "--penalty", metavar="P", type=int, nargs=None,
                        help="Penalty to administer when opening or expanding gap. Default is 1.")

    parser.add_argument("-s1", "--sequence1", metavar="S", type=str,
                        nargs=None, help="Sequence number one.")
    parser.add_argument("-s2", "--sequence2", metavar="S", type=str,
                        nargs=None, help="Sequence number two.")

    parser.add_argument("-f1", "--path1", metavar="F", type=str,
                        nargs=None, help="Path to file containing Sequence number one.")
    parser.add_argument("-f2", "--path2", metavar="F", type=str,
                        nargs=None, help="Path to file containing Sequence number two.")

    args = parser.parse_args()

    seq1 = args.sequence1
    seq2 = args.sequence2
    p = args.penalty

    if not args.sequence1:
        seq1 = read_file_content(args.path1)

    if not args.sequence2:
        seq2 = read_file_content(args.path2)

    return p, seq1, seq2


if __name__ == "__main__":
    # penalty, seq1, seq2 = parse_args()
    seq1 = "TECTEA"
    seq2 = "CCTEC"
    penalty = 1
    sw = SmithWaterman(penalty if penalty else 1, seq1, seq2)

    sw.solve()

    traces = sw.traceback()
    print(traces)
