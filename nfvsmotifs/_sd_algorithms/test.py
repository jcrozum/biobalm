from expand_source_SCCs import *


bn = BooleanNetwork.from_file("nfvsmotifs/_sd_algorithms/test.bnet")

source_nodes = find_source_nodes(bn)

print(source_nodes)

perc_space, _ = percolate_space(bn, {}, strict_percolation=False)
perc_bn = percolate_network(bn, perc_space)

source_nodes = find_source_nodes(perc_bn)

print(source_nodes)

bin_values_list = it.product(range(2), repeat=len(source_nodes)) 
for bin_values in bin_values_list:
    source_comb = dict(zip(source_nodes, bin_values))
    print(source_comb)

    sub_space = source_comb
    sub_space.update(perc_space)

    print(sub_space)