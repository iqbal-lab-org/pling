from unittest import TestCase
from unittest.mock import patch
from pathlib import Path
import ast
import itertools
import math
import tempfile
import numpy as np
from pling.batching.get_batches import get_labels, append_pair, get_pairs, main


def read_batches(output_dir):
    """Returns the pairs in each batch file, ordered by batch number."""
    batch_files = sorted(Path(output_dir).glob("batch_*.txt"), key=lambda p: int(p.stem.split("_")[1]))
    return [[tuple(ast.literal_eval(line)) for line in f.read_text().splitlines()] for f in batch_files]

def as_pair_set(pairs):
    return {frozenset(pair) for pair in pairs}

def expected_pairs(genomes, prev_genomes=[]):
    """All unordered pairs, except those already compared within the same previous typing."""
    return {frozenset(pair) for pair in itertools.combinations(genomes, 2)
            if not any(set(pair) <= set(typing) for typing in prev_genomes)}

class Test_get_pairs(TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.output_dir = Path(self.tmp.name)/"batches"
        self.genomes = [f"g{i}" for i in range(7)]

    def tearDown(self):
        self.tmp.cleanup()

    def test_no_previous_typing_gives_all_pairs(self):
        n_pairs = get_pairs(self.genomes, 4, self.output_dir, None, False, False)
        batches = read_batches(self.output_dir)
        pairs = list(itertools.chain(*batches))
        self.assertEqual(n_pairs, 21)
        self.assertEqual(len(pairs), 21)
        self.assertEqual(as_pair_set(pairs), expected_pairs(self.genomes))
        self.assertTrue(all(genome_1 != genome_2 for genome_1, genome_2 in pairs))

    def test_batch_sizes(self):
        get_pairs(self.genomes, 4, self.output_dir, None, False, False)
        batches = read_batches(self.output_dir)
        self.assertEqual(len(batches), math.ceil(21/4))
        self.assertEqual([len(batch) for batch in batches], [4, 4, 4, 4, 4, 1])

    def test_batch_size_divides_number_of_pairs(self):
        get_pairs(self.genomes, 7, self.output_dir, None, False, False)
        self.assertEqual([len(batch) for batch in read_batches(self.output_dir)], [7, 7, 7])

    def test_previous_typings(self):
        prev_genomes = [["g0", "g1", "g2"], ["g3", "g4"]]
        new_genomes = ["g5", "g6"]
        n_pairs = get_pairs(self.genomes, 4, self.output_dir, None, False, False, prev_genomes=prev_genomes, new_genomes=new_genomes)
        pairs = list(itertools.chain(*read_batches(self.output_dir)))
        self.assertEqual(n_pairs, 17)
        self.assertEqual(len(pairs), 17)
        self.assertEqual(as_pair_set(pairs), expected_pairs(self.genomes, prev_genomes))

    def test_previous_typings_without_new_genomes(self):
        prev_genomes = [["g0", "g1", "g2"], ["g3", "g4", "g5", "g6"]]
        n_pairs = get_pairs(self.genomes, 4, self.output_dir, None, False, False, prev_genomes=prev_genomes, new_genomes=[])
        pairs = list(itertools.chain(*read_batches(self.output_dir)))
        self.assertEqual(n_pairs, 12)
        self.assertEqual(as_pair_set(pairs), expected_pairs(self.genomes, prev_genomes))

    def test_single_previous_typing_with_new_genomes(self):
        prev_genomes = [["g0", "g1", "g2", "g3", "g4"]]
        new_genomes = ["g5", "g6"]
        n_pairs = get_pairs(self.genomes, 4, self.output_dir, None, False, False, prev_genomes=prev_genomes, new_genomes=new_genomes)
        pairs = list(itertools.chain(*read_batches(self.output_dir)))
        self.assertEqual(n_pairs, 11)
        self.assertEqual(as_pair_set(pairs), expected_pairs(self.genomes, prev_genomes))

    def test_no_pairs_left(self):
        prev_genomes = [self.genomes]
        n_pairs = get_pairs(self.genomes, 4, self.output_dir, None, False, False, prev_genomes=prev_genomes, new_genomes=[])
        self.assertEqual(n_pairs, 0)
        self.assertEqual(read_batches(self.output_dir), [])

    def test_sourmash_filters_pairs_and_writes_containment(self):
        genomes = ["a", "b", "c"]
        # containment similarity, distance is 1-similarity
        smash_matrix = np.array([[1, 0.9, 0.1],
                                 [0.9, 1, 0.8],
                                 [0.1, 0.8, 1]])
        containmentpath = Path(self.tmp.name)/"containment"/"distances.tsv"
        n_pairs = get_pairs(genomes, 4, self.output_dir, containmentpath, True, False, smash_matrix, 0.5)
        pairs = list(itertools.chain(*read_batches(self.output_dir)))
        self.assertEqual(n_pairs, 2)
        self.assertEqual(as_pair_set(pairs), {frozenset(["a", "b"]), frozenset(["b", "c"])})
        # only pairs filtered out by sourmash get their containment distance written
        lines = containmentpath.read_text().splitlines()
        self.assertEqual(lines[0], "plasmid_1\tplasmid_2\tdistance")
        self.assertEqual(len(lines), 2)
        plasmid_1, plasmid_2, dist = lines[1].split("\t")
        self.assertEqual({plasmid_1, plasmid_2}, {"a", "c"})
        self.assertAlmostEqual(float(dist), 0.9)

    def test_sourmash_only_writes_containment_for_batched_pairs(self):
        genomes = ["a", "b", "c"]
        smash_matrix = np.array([[1, 0.9, 0.1],
                                 [0.9, 1, 0.8],
                                 [0.1, 0.8, 1]])
        containmentpath = Path(self.tmp.name)/"distances.tsv"
        n_pairs = get_pairs(genomes, 4, self.output_dir, containmentpath, False, True, smash_matrix, 0.5)
        self.assertEqual(n_pairs, 3)
        lines = containmentpath.read_text().splitlines()[1:]
        self.assertEqual({frozenset(line.split("\t")[:2]) for line in lines}, expected_pairs(genomes))


class Test_main(TestCase):
    def test_batching_info_without_sourmash(self):
        with tempfile.TemporaryDirectory() as tmp:
            genomes_list = Path(tmp)/"genomes.txt"
            genomes_list.write_text("".join(f"data/g{i}.fna\n" for i in range(5)))
            outputpath = Path(tmp)/"out"
            argv = ["get_batches.py", "--genomes_list", str(genomes_list), "--batch_size", "3", "--outputpath", str(outputpath)]
            with patch("sys.argv", argv):
                main()
            batch_size, number_of_batches = (outputpath/"batches"/"batching_info.txt").read_text().split("\n")
            self.assertEqual(int(batch_size), 3)
            self.assertEqual(int(number_of_batches), math.ceil(10/3))
            self.assertEqual(len(read_batches(outputpath/"batches")), int(number_of_batches))
