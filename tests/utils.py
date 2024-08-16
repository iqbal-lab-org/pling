import unittest

# Instance of the TestCase to use the assertEqual method
tc = unittest.TestCase()


def read_file(path):
    """Helper method to read a file."""
    with open(path, 'rb') as f:
        return f.read()


def assert_files_are_identical(path_to_first_file, path_to_second_file):
    first_file_content = read_file(path_to_first_file)
    second_file_content = read_file(path_to_second_file)
    tc.assertEqual(first_file_content, second_file_content)

def get_num_communities(filepath):
    with open(filepath) as communities_fh:
        num = len(communities_fh.readlines())
    return num

def get_batch_num(filepath):
    with open(filepath) as f:
        next(f)
        num = int(f.readline().strip())
    return num

def assert_containment(test, test_dir, integerisation):
    assert_files_are_identical(f"{test_dir}/{test}/out/{integerisation}/containment/all_pairs_containment_distance.tsv",
                               f"{test_dir}/{test}/truth/{integerisation}/containment/all_pairs_containment_distance.tsv")
    assert_files_are_identical(f"{test_dir}/{test}/out/{integerisation}/containment/containment_communities/objects/communities.tsv",
                               f"{test_dir}/{test}/truth/{integerisation}/containment/containment_communities/objects/communities.tsv")
    assert_files_are_identical(f"{test_dir}/{test}/out/{integerisation}/containment/containment_communities/objects/communities.txt",
                               f"{test_dir}/{test}/truth/{integerisation}/containment/containment_communities/objects/communities.txt")

def assert_align_unimogs(test, test_dir, batch_num):
    for i in range(batch_num):
        assert_files_are_identical(f"{test_dir}/{test}/out/align/unimogs/batch_{i}_align.unimog",
                                   f"{test_dir}/{test}/truth/align/unimogs/batch_{i}_align.unimog")
        assert_files_are_identical(f"{test_dir}/{test}/out/align/unimogs/batch_{i}_map.txt",
                                   f"{test_dir}/{test}/truth/align/unimogs/batch_{i}_map.txt")

def assert_end_to_end(test, dir):
    #containment communities
    assert_containment(test, dir, "align")
    #unimogs
    batch_num = get_batch_num(f"{dir}/{test}/out/align/batches/batching_info.txt")
    assert_align_unimogs(test, dir, batch_num)
    #DCJ distance matrix
    assert_files_are_identical(f"{dir}/{test}/out/align/all_plasmids_distances.tsv",
                               f"{dir}/{test}/truth/align/all_plasmids_distances.tsv")
    #DCJ communities
    assert_files_are_identical(f"{dir}/{test}/out/align/dcj_thresh_4_graph/objects/typing.tsv",
                               f"{dir}/{test}/truth/align/dcj_thresh_4_graph/objects/typing.tsv")
    assert_files_are_identical(f"{dir}/{test}/out/align/dcj_thresh_4_graph/objects/hub_plasmids.csv",
                               f"{dir}/{test}/truth/align/dcj_thresh_4_graph/objects/hub_plasmids.csv")
