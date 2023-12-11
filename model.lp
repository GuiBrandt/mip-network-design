\ Model Um Problema de Network Design
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  3767 circuit_arc[4,5] + 3767 circuit_arc[5,4] + 2760 circuit_arc[4,6]
   + 2760 circuit_arc[6,4] + 4411 circuit_arc[4,7] + 4411 circuit_arc[7,4]
   + 2844 circuit_arc[4,8] + 2844 circuit_arc[8,4] + 5088 circuit_arc[4,9]
   + 5088 circuit_arc[9,4] + 3187 circuit_arc[3,9] + 3187 circuit_arc[9,3]
   + 6066 circuit_arc[3,8] + 6066 circuit_arc[8,3] + 2046 circuit_arc[3,7]
   + 2046 circuit_arc[7,3] + 6406 circuit_arc[3,6] + 6406 circuit_arc[6,3]
   + 6644 circuit_arc[3,5] + 6644 circuit_arc[5,3] + 3868 circuit_arc[3,4]
   + 3868 circuit_arc[4,3] + 2727 circuit_arc[5,6] + 2727 circuit_arc[6,5]
   + 7793 circuit_arc[5,7] + 7793 circuit_arc[7,5] + 1485 circuit_arc[5,8]
   + 1485 circuit_arc[8,5] + 6066 circuit_arc[5,9] + 6066 circuit_arc[9,5]
   + 3651 circuit_arc[2,9] + 3651 circuit_arc[9,2] + 2726 circuit_arc[2,8]
   + 2726 circuit_arc[8,2] + 4649 circuit_arc[2,7] + 4649 circuit_arc[7,2]
   + 3288 circuit_arc[2,6] + 3288 circuit_arc[6,2] + 3600 circuit_arc[2,5]
   + 3600 circuit_arc[5,2] + 1724 circuit_arc[2,4] + 1724 circuit_arc[4,2]
   + 3435 circuit_arc[2,3] + 3435 circuit_arc[3,2] + 7382 circuit_arc[6,7]
   + 7382 circuit_arc[7,6] + 1504 circuit_arc[6,8] + 1504 circuit_arc[8,6]
   + 7070 circuit_arc[6,9] + 7070 circuit_arc[9,6] + 5408 circuit_arc[1,9]
   + 5408 circuit_arc[9,1] + 4410 circuit_arc[1,8] + 4410 circuit_arc[8,1]
   + 3652 circuit_arc[1,7] + 3652 circuit_arc[7,1] + 4117 circuit_arc[1,6]
   + 4117 circuit_arc[6,1] + 5428 circuit_arc[1,5] + 5428 circuit_arc[5,1]
   + 1478 circuit_arc[1,4] + 1478 circuit_arc[4,1] + 3111 circuit_arc[1,3]
   + 3111 circuit_arc[3,1] + 2677 circuit_arc[1,2] + 2677 circuit_arc[2,1]
   + 7350 circuit_arc[7,8] + 7350 circuit_arc[8,7] + 4378 circuit_arc[7,9]
   + 4378 circuit_arc[9,7] + 5555 circuit_arc[0,9] + 5555 circuit_arc[9,0]
   + 1693 circuit_arc[0,8] + 1693 circuit_arc[8,0] + 6110 circuit_arc[0,7]
   + 6110 circuit_arc[7,0] + 1481 circuit_arc[0,6] + 1481 circuit_arc[6,0]
   + 2473 circuit_arc[0,5] + 2473 circuit_arc[5,0] + 1889 circuit_arc[0,4]
   + 1889 circuit_arc[4,0] + 4952 circuit_arc[0,3] + 4952 circuit_arc[3,0]
   + 2246 circuit_arc[0,2] + 2246 circuit_arc[2,0] + 3066 circuit_arc[0,1]
   + 3066 circuit_arc[1,0] + 5821 circuit_arc[8,9] + 5821 circuit_arc[9,8]
   + 3767 star_edge[4,5] + 2760 star_edge[4,6] + 4411 star_edge[4,7]
   + 2844 star_edge[4,8] + 5088 star_edge[4,9] + 3187 star_edge[3,9]
   + 6066 star_edge[3,8] + 2046 star_edge[3,7] + 6406 star_edge[3,6]
   + 6644 star_edge[3,5] + 3868 star_edge[3,4] + 2727 star_edge[5,6]
   + 7793 star_edge[5,7] + 1485 star_edge[5,8] + 6066 star_edge[5,9]
   + 3651 star_edge[2,9] + 2726 star_edge[2,8] + 4649 star_edge[2,7]
   + 3288 star_edge[2,6] + 3600 star_edge[2,5] + 1724 star_edge[2,4]
   + 3435 star_edge[2,3] + 7382 star_edge[6,7] + 1504 star_edge[6,8]
   + 7070 star_edge[6,9] + 5408 star_edge[1,9] + 4410 star_edge[1,8]
   + 3652 star_edge[1,7] + 4117 star_edge[1,6] + 5428 star_edge[1,5]
   + 1478 star_edge[1,4] + 3111 star_edge[1,3] + 2677 star_edge[1,2]
   + 7350 star_edge[7,8] + 4378 star_edge[7,9] + 5555 star_edge[0,9]
   + 1693 star_edge[0,8] + 6110 star_edge[0,7] + 1481 star_edge[0,6]
   + 2473 star_edge[0,5] + 1889 star_edge[0,4] + 4952 star_edge[0,3]
   + 2246 star_edge[0,2] + 3066 star_edge[0,1] + 5821 star_edge[8,9]
Subject To
 Partition_0_-_Used_with_node_9: - partition[9,0] + partition_used[0] >= 0
 Partition_0_-_Used_with_node_8: - partition[8,0] + partition_used[0] >= 0
 Partition_0_-_Used_with_node_7: - partition[7,0] + partition_used[0] >= 0
 Partition_0_-_Used_with_node_6: - partition[6,0] + partition_used[0] >= 0
 Partition_0_-_Used_with_node_5: - partition[5,0] + partition_used[0] >= 0
 Partition_0_-_Used_with_node_4: - partition[4,0] + partition_used[0] >= 0
 Partition_0_-_Used_with_node_3: - partition[3,0] + partition_used[0] >= 0
 Partition_0_-_Used_with_node_2: - partition[2,0] + partition_used[0] >= 0
 Partition_0_-_Used_with_node_1: - partition[1,0] + partition_used[0] >= 0
 Partition_0_-_Used_with_node_0: - partition[0,0] + partition_used[0] >= 0
 Partition_1_-_Used_with_node_9: - partition[9,1] + partition_used[1] >= 0
 Partition_1_-_Used_with_node_8: - partition[8,1] + partition_used[1] >= 0
 Partition_1_-_Used_with_node_7: - partition[7,1] + partition_used[1] >= 0
 Partition_1_-_Used_with_node_6: - partition[6,1] + partition_used[1] >= 0
 Partition_1_-_Used_with_node_5: - partition[5,1] + partition_used[1] >= 0
 Partition_1_-_Used_with_node_4: - partition[4,1] + partition_used[1] >= 0
 Partition_1_-_Used_with_node_3: - partition[3,1] + partition_used[1] >= 0
 Partition_1_-_Used_with_node_2: - partition[2,1] + partition_used[1] >= 0
 Partition_1_-_Used_with_node_1: - partition[1,1] + partition_used[1] >= 0
 Partition_1_-_Used_with_node_0: - partition[0,1] + partition_used[1] >= 0
 Partition_2_-_Used_with_node_9: - partition[9,2] + partition_used[2] >= 0
 Partition_2_-_Used_with_node_8: - partition[8,2] + partition_used[2] >= 0
 Partition_2_-_Used_with_node_7: - partition[7,2] + partition_used[2] >= 0
 Partition_2_-_Used_with_node_6: - partition[6,2] + partition_used[2] >= 0
 Partition_2_-_Used_with_node_5: - partition[5,2] + partition_used[2] >= 0
 Partition_2_-_Used_with_node_4: - partition[4,2] + partition_used[2] >= 0
 Partition_2_-_Used_with_node_3: - partition[3,2] + partition_used[2] >= 0
 Partition_2_-_Used_with_node_2: - partition[2,2] + partition_used[2] >= 0
 Partition_2_-_Used_with_node_1: - partition[1,2] + partition_used[2] >= 0
 Partition_2_-_Used_with_node_0: - partition[0,2] + partition_used[2] >= 0
 Partition_3_-_Used_with_node_9: - partition[9,3] + partition_used[3] >= 0
 Partition_3_-_Used_with_node_8: - partition[8,3] + partition_used[3] >= 0
 Partition_3_-_Used_with_node_7: - partition[7,3] + partition_used[3] >= 0
 Partition_3_-_Used_with_node_6: - partition[6,3] + partition_used[3] >= 0
 Partition_3_-_Used_with_node_5: - partition[5,3] + partition_used[3] >= 0
 Partition_3_-_Used_with_node_4: - partition[4,3] + partition_used[3] >= 0
 Partition_3_-_Used_with_node_3: - partition[3,3] + partition_used[3] >= 0
 Partition_3_-_Used_with_node_2: - partition[2,3] + partition_used[3] >= 0
 Partition_3_-_Used_with_node_1: - partition[1,3] + partition_used[3] >= 0
 Partition_3_-_Used_with_node_0: - partition[0,3] + partition_used[3] >= 0
 Partition_4_-_Used_with_node_9: - partition[9,4] + partition_used[4] >= 0
 Partition_4_-_Used_with_node_8: - partition[8,4] + partition_used[4] >= 0
 Partition_4_-_Used_with_node_7: - partition[7,4] + partition_used[4] >= 0
 Partition_4_-_Used_with_node_6: - partition[6,4] + partition_used[4] >= 0
 Partition_4_-_Used_with_node_5: - partition[5,4] + partition_used[4] >= 0
 Partition_4_-_Used_with_node_4: - partition[4,4] + partition_used[4] >= 0
 Partition_4_-_Used_with_node_3: - partition[3,4] + partition_used[4] >= 0
 Partition_4_-_Used_with_node_2: - partition[2,4] + partition_used[4] >= 0
 Partition_4_-_Used_with_node_1: - partition[1,4] + partition_used[4] >= 0
 Partition_4_-_Used_with_node_0: - partition[0,4] + partition_used[4] >= 0
 Partition_5_-_Used_with_node_9: - partition[9,5] + partition_used[5] >= 0
 Partition_5_-_Used_with_node_8: - partition[8,5] + partition_used[5] >= 0
 Partition_5_-_Used_with_node_7: - partition[7,5] + partition_used[5] >= 0
 Partition_5_-_Used_with_node_6: - partition[6,5] + partition_used[5] >= 0
 Partition_5_-_Used_with_node_5: - partition[5,5] + partition_used[5] >= 0
 Partition_5_-_Used_with_node_4: - partition[4,5] + partition_used[5] >= 0
 Partition_5_-_Used_with_node_3: - partition[3,5] + partition_used[5] >= 0
 Partition_5_-_Used_with_node_2: - partition[2,5] + partition_used[5] >= 0
 Partition_5_-_Used_with_node_1: - partition[1,5] + partition_used[5] >= 0
 Partition_5_-_Used_with_node_0: - partition[0,5] + partition_used[5] >= 0
 Partition_6_-_Used_with_node_9: - partition[9,6] + partition_used[6] >= 0
 Partition_6_-_Used_with_node_8: - partition[8,6] + partition_used[6] >= 0
 Partition_6_-_Used_with_node_7: - partition[7,6] + partition_used[6] >= 0
 Partition_6_-_Used_with_node_6: - partition[6,6] + partition_used[6] >= 0
 Partition_6_-_Used_with_node_5: - partition[5,6] + partition_used[6] >= 0
 Partition_6_-_Used_with_node_4: - partition[4,6] + partition_used[6] >= 0
 Partition_6_-_Used_with_node_3: - partition[3,6] + partition_used[6] >= 0
 Partition_6_-_Used_with_node_2: - partition[2,6] + partition_used[6] >= 0
 Partition_6_-_Used_with_node_1: - partition[1,6] + partition_used[6] >= 0
 Partition_6_-_Used_with_node_0: - partition[0,6] + partition_used[6] >= 0
 Partition_7_-_Used_with_node_9: - partition[9,7] + partition_used[7] >= 0
 Partition_7_-_Used_with_node_8: - partition[8,7] + partition_used[7] >= 0
 Partition_7_-_Used_with_node_7: - partition[7,7] + partition_used[7] >= 0
 Partition_7_-_Used_with_node_6: - partition[6,7] + partition_used[7] >= 0
 Partition_7_-_Used_with_node_5: - partition[5,7] + partition_used[7] >= 0
 Partition_7_-_Used_with_node_4: - partition[4,7] + partition_used[7] >= 0
 Partition_7_-_Used_with_node_3: - partition[3,7] + partition_used[7] >= 0
 Partition_7_-_Used_with_node_2: - partition[2,7] + partition_used[7] >= 0
 Partition_7_-_Used_with_node_1: - partition[1,7] + partition_used[7] >= 0
 Partition_7_-_Used_with_node_0: - partition[0,7] + partition_used[7] >= 0
 Partition_8_-_Used_with_node_9: - partition[9,8] + partition_used[8] >= 0
 Partition_8_-_Used_with_node_8: - partition[8,8] + partition_used[8] >= 0
 Partition_8_-_Used_with_node_7: - partition[7,8] + partition_used[8] >= 0
 Partition_8_-_Used_with_node_6: - partition[6,8] + partition_used[8] >= 0
 Partition_8_-_Used_with_node_5: - partition[5,8] + partition_used[8] >= 0
 Partition_8_-_Used_with_node_4: - partition[4,8] + partition_used[8] >= 0
 Partition_8_-_Used_with_node_3: - partition[3,8] + partition_used[8] >= 0
 Partition_8_-_Used_with_node_2: - partition[2,8] + partition_used[8] >= 0
 Partition_8_-_Used_with_node_1: - partition[1,8] + partition_used[8] >= 0
 Partition_8_-_Used_with_node_0: - partition[0,8] + partition_used[8] >= 0
 Partition_9_-_Used_with_node_9: - partition[9,9] + partition_used[9] >= 0
 Partition_9_-_Used_with_node_8: - partition[8,9] + partition_used[9] >= 0
 Partition_9_-_Used_with_node_7: - partition[7,9] + partition_used[9] >= 0
 Partition_9_-_Used_with_node_6: - partition[6,9] + partition_used[9] >= 0
 Partition_9_-_Used_with_node_5: - partition[5,9] + partition_used[9] >= 0
 Partition_9_-_Used_with_node_4: - partition[4,9] + partition_used[9] >= 0
 Partition_9_-_Used_with_node_3: - partition[3,9] + partition_used[9] >= 0
 Partition_9_-_Used_with_node_2: - partition[2,9] + partition_used[9] >= 0
 Partition_9_-_Used_with_node_1: - partition[1,9] + partition_used[9] >= 0
 Partition_9_-_Used_with_node_0: - partition[0,9] + partition_used[9] >= 0
 At_least_3_partitions: partition_used[0] + partition_used[1]
   + partition_used[2] + partition_used[3] + partition_used[4]
   + partition_used[5] + partition_used[6] + partition_used[7]
   + partition_used[8] + partition_used[9] >= 3
 Partition_0_-_Capacity_limit: 46 partition[9,0] + 68 partition[8,0]
   + 26 partition[7,0] + 56 partition[6,0] + 78 partition[5,0]
   + 29 partition[4,0] + 54 partition[3,0] + 31 partition[2,0]
   + 51 partition[1,0] + 73 partition[0,0] <= 100
 Partition_1_-_Capacity_limit: 46 partition[9,1] + 68 partition[8,1]
   + 26 partition[7,1] + 56 partition[6,1] + 78 partition[5,1]
   + 29 partition[4,1] + 54 partition[3,1] + 31 partition[2,1]
   + 51 partition[1,1] + 73 partition[0,1] <= 100
 Partition_2_-_Capacity_limit: 46 partition[9,2] + 68 partition[8,2]
   + 26 partition[7,2] + 56 partition[6,2] + 78 partition[5,2]
   + 29 partition[4,2] + 54 partition[3,2] + 31 partition[2,2]
   + 51 partition[1,2] + 73 partition[0,2] <= 100
 Partition_3_-_Capacity_limit: 46 partition[9,3] + 68 partition[8,3]
   + 26 partition[7,3] + 56 partition[6,3] + 78 partition[5,3]
   + 29 partition[4,3] + 54 partition[3,3] + 31 partition[2,3]
   + 51 partition[1,3] + 73 partition[0,3] <= 100
 Partition_4_-_Capacity_limit: 46 partition[9,4] + 68 partition[8,4]
   + 26 partition[7,4] + 56 partition[6,4] + 78 partition[5,4]
   + 29 partition[4,4] + 54 partition[3,4] + 31 partition[2,4]
   + 51 partition[1,4] + 73 partition[0,4] <= 100
 Partition_5_-_Capacity_limit: 46 partition[9,5] + 68 partition[8,5]
   + 26 partition[7,5] + 56 partition[6,5] + 78 partition[5,5]
   + 29 partition[4,5] + 54 partition[3,5] + 31 partition[2,5]
   + 51 partition[1,5] + 73 partition[0,5] <= 100
 Partition_6_-_Capacity_limit: 46 partition[9,6] + 68 partition[8,6]
   + 26 partition[7,6] + 56 partition[6,6] + 78 partition[5,6]
   + 29 partition[4,6] + 54 partition[3,6] + 31 partition[2,6]
   + 51 partition[1,6] + 73 partition[0,6] <= 100
 Partition_7_-_Capacity_limit: 46 partition[9,7] + 68 partition[8,7]
   + 26 partition[7,7] + 56 partition[6,7] + 78 partition[5,7]
   + 29 partition[4,7] + 54 partition[3,7] + 31 partition[2,7]
   + 51 partition[1,7] + 73 partition[0,7] <= 100
 Partition_8_-_Capacity_limit: 46 partition[9,8] + 68 partition[8,8]
   + 26 partition[7,8] + 56 partition[6,8] + 78 partition[5,8]
   + 29 partition[4,8] + 54 partition[3,8] + 31 partition[2,8]
   + 51 partition[1,8] + 73 partition[0,8] <= 100
 Partition_9_-_Capacity_limit: 46 partition[9,9] + 68 partition[8,9]
   + 26 partition[7,9] + 56 partition[6,9] + 78 partition[5,9]
   + 29 partition[4,9] + 54 partition[3,9] + 31 partition[2,9]
   + 51 partition[1,9] + 73 partition[0,9] <= 100
 Node_9_-_Exactly_one_partition: partition[9,0] + partition[9,1]
   + partition[9,2] + partition[9,3] + partition[9,4] + partition[9,5]
   + partition[9,6] + partition[9,7] + partition[9,8] + partition[9,9] = 1
 Node_8_-_Exactly_one_partition: partition[8,0] + partition[8,1]
   + partition[8,2] + partition[8,3] + partition[8,4] + partition[8,5]
   + partition[8,6] + partition[8,7] + partition[8,8] + partition[8,9] = 1
 Node_7_-_Exactly_one_partition: partition[7,0] + partition[7,1]
   + partition[7,2] + partition[7,3] + partition[7,4] + partition[7,5]
   + partition[7,6] + partition[7,7] + partition[7,8] + partition[7,9] = 1
 Node_6_-_Exactly_one_partition: partition[6,0] + partition[6,1]
   + partition[6,2] + partition[6,3] + partition[6,4] + partition[6,5]
   + partition[6,6] + partition[6,7] + partition[6,8] + partition[6,9] = 1
 Node_5_-_Exactly_one_partition: partition[5,0] + partition[5,1]
   + partition[5,2] + partition[5,3] + partition[5,4] + partition[5,5]
   + partition[5,6] + partition[5,7] + partition[5,8] + partition[5,9] = 1
 Node_4_-_Exactly_one_partition: partition[4,0] + partition[4,1]
   + partition[4,2] + partition[4,3] + partition[4,4] + partition[4,5]
   + partition[4,6] + partition[4,7] + partition[4,8] + partition[4,9] = 1
 Node_3_-_Exactly_one_partition: partition[3,0] + partition[3,1]
   + partition[3,2] + partition[3,3] + partition[3,4] + partition[3,5]
   + partition[3,6] + partition[3,7] + partition[3,8] + partition[3,9] = 1
 Node_2_-_Exactly_one_partition: partition[2,0] + partition[2,1]
   + partition[2,2] + partition[2,3] + partition[2,4] + partition[2,5]
   + partition[2,6] + partition[2,7] + partition[2,8] + partition[2,9] = 1
 Node_1_-_Exactly_one_partition: partition[1,0] + partition[1,1]
   + partition[1,2] + partition[1,3] + partition[1,4] + partition[1,5]
   + partition[1,6] + partition[1,7] + partition[1,8] + partition[1,9] = 1
 Node_0_-_Exactly_one_partition: partition[0,0] + partition[0,1]
   + partition[0,2] + partition[0,3] + partition[0,4] + partition[0,5]
   + partition[0,6] + partition[0,7] + partition[0,8] + partition[0,9] = 1
 Partition_1_-_Symmetry_breaking: - partition_used[0] + partition_used[1]
   <= 0
 Partition_2_-_Symmetry_breaking: - partition_used[1] + partition_used[2]
   <= 0
 Partition_3_-_Symmetry_breaking: - partition_used[2] + partition_used[3]
   <= 0
 Partition_4_-_Symmetry_breaking: - partition_used[3] + partition_used[4]
   <= 0
 Partition_5_-_Symmetry_breaking: - partition_used[4] + partition_used[5]
   <= 0
 Partition_6_-_Symmetry_breaking: - partition_used[5] + partition_used[6]
   <= 0
 Partition_7_-_Symmetry_breaking: - partition_used[6] + partition_used[7]
   <= 0
 Partition_8_-_Symmetry_breaking: - partition_used[7] + partition_used[8]
   <= 0
 Partition_9_-_Symmetry_breaking: - partition_used[8] + partition_used[9]
   <= 0
 Partition_0_-_Star_edge_{4,_5}: - 2 partition[5,0]
   - circuit_partition_node[4,0] + star_edge[4,5] >= -2
 Partition_0_-_Star_edge_{4,_5}: - 2 partition[4,0]
   - circuit_partition_node[5,0] + star_edge[4,5] >= -2
 Partition_0_-_Star_edge_{4,_6}: - 2 partition[6,0]
   - circuit_partition_node[4,0] + star_edge[4,6] >= -2
 Partition_0_-_Star_edge_{4,_6}: - 2 partition[4,0]
   - circuit_partition_node[6,0] + star_edge[4,6] >= -2
 Partition_0_-_Star_edge_{4,_7}: - 2 partition[7,0]
   - circuit_partition_node[4,0] + star_edge[4,7] >= -2
 Partition_0_-_Star_edge_{4,_7}: - 2 partition[4,0]
   - circuit_partition_node[7,0] + star_edge[4,7] >= -2
 Partition_0_-_Star_edge_{4,_8}: - 2 partition[8,0]
   - circuit_partition_node[4,0] + star_edge[4,8] >= -2
 Partition_0_-_Star_edge_{4,_8}: - 2 partition[4,0]
   - circuit_partition_node[8,0] + star_edge[4,8] >= -2
 Partition_0_-_Star_edge_{4,_9}: - 2 partition[9,0]
   - circuit_partition_node[4,0] + star_edge[4,9] >= -2
 Partition_0_-_Star_edge_{4,_9}: - 2 partition[4,0]
   - circuit_partition_node[9,0] + star_edge[4,9] >= -2
 Partition_0_-_Star_edge_{3,_9}: - 2 partition[9,0]
   - circuit_partition_node[3,0] + star_edge[3,9] >= -2
 Partition_0_-_Star_edge_{3,_9}: - 2 partition[3,0]
   - circuit_partition_node[9,0] + star_edge[3,9] >= -2
 Partition_0_-_Star_edge_{3,_8}: - 2 partition[8,0]
   - circuit_partition_node[3,0] + star_edge[3,8] >= -2
 Partition_0_-_Star_edge_{3,_8}: - 2 partition[3,0]
   - circuit_partition_node[8,0] + star_edge[3,8] >= -2
 Partition_0_-_Star_edge_{3,_7}: - 2 partition[7,0]
   - circuit_partition_node[3,0] + star_edge[3,7] >= -2
 Partition_0_-_Star_edge_{3,_7}: - 2 partition[3,0]
   - circuit_partition_node[7,0] + star_edge[3,7] >= -2
 Partition_0_-_Star_edge_{3,_6}: - 2 partition[6,0]
   - circuit_partition_node[3,0] + star_edge[3,6] >= -2
 Partition_0_-_Star_edge_{3,_6}: - 2 partition[3,0]
   - circuit_partition_node[6,0] + star_edge[3,6] >= -2
 Partition_0_-_Star_edge_{3,_5}: - 2 partition[5,0]
   - circuit_partition_node[3,0] + star_edge[3,5] >= -2
 Partition_0_-_Star_edge_{3,_5}: - 2 partition[3,0]
   - circuit_partition_node[5,0] + star_edge[3,5] >= -2
 Partition_0_-_Star_edge_{3,_4}: - 2 partition[4,0]
   - circuit_partition_node[3,0] + star_edge[3,4] >= -2
 Partition_0_-_Star_edge_{3,_4}: - 2 partition[3,0]
   - circuit_partition_node[4,0] + star_edge[3,4] >= -2
 Partition_0_-_Star_edge_{5,_6}: - 2 partition[6,0]
   - circuit_partition_node[5,0] + star_edge[5,6] >= -2
 Partition_0_-_Star_edge_{5,_6}: - 2 partition[5,0]
   - circuit_partition_node[6,0] + star_edge[5,6] >= -2
 Partition_0_-_Star_edge_{5,_7}: - 2 partition[7,0]
   - circuit_partition_node[5,0] + star_edge[5,7] >= -2
 Partition_0_-_Star_edge_{5,_7}: - 2 partition[5,0]
   - circuit_partition_node[7,0] + star_edge[5,7] >= -2
 Partition_0_-_Star_edge_{5,_8}: - 2 partition[8,0]
   - circuit_partition_node[5,0] + star_edge[5,8] >= -2
 Partition_0_-_Star_edge_{5,_8}: - 2 partition[5,0]
   - circuit_partition_node[8,0] + star_edge[5,8] >= -2
 Partition_0_-_Star_edge_{5,_9}: - 2 partition[9,0]
   - circuit_partition_node[5,0] + star_edge[5,9] >= -2
 Partition_0_-_Star_edge_{5,_9}: - 2 partition[5,0]
   - circuit_partition_node[9,0] + star_edge[5,9] >= -2
 Partition_0_-_Star_edge_{2,_9}: - 2 partition[9,0]
   - circuit_partition_node[2,0] + star_edge[2,9] >= -2
 Partition_0_-_Star_edge_{2,_9}: - 2 partition[2,0]
   - circuit_partition_node[9,0] + star_edge[2,9] >= -2
 Partition_0_-_Star_edge_{2,_8}: - 2 partition[8,0]
   - circuit_partition_node[2,0] + star_edge[2,8] >= -2
 Partition_0_-_Star_edge_{2,_8}: - 2 partition[2,0]
   - circuit_partition_node[8,0] + star_edge[2,8] >= -2
 Partition_0_-_Star_edge_{2,_7}: - 2 partition[7,0]
   - circuit_partition_node[2,0] + star_edge[2,7] >= -2
 Partition_0_-_Star_edge_{2,_7}: - 2 partition[2,0]
   - circuit_partition_node[7,0] + star_edge[2,7] >= -2
 Partition_0_-_Star_edge_{2,_6}: - 2 partition[6,0]
   - circuit_partition_node[2,0] + star_edge[2,6] >= -2
 Partition_0_-_Star_edge_{2,_6}: - 2 partition[2,0]
   - circuit_partition_node[6,0] + star_edge[2,6] >= -2
 Partition_0_-_Star_edge_{2,_5}: - 2 partition[5,0]
   - circuit_partition_node[2,0] + star_edge[2,5] >= -2
 Partition_0_-_Star_edge_{2,_5}: - 2 partition[2,0]
   - circuit_partition_node[5,0] + star_edge[2,5] >= -2
 Partition_0_-_Star_edge_{2,_4}: - 2 partition[4,0]
   - circuit_partition_node[2,0] + star_edge[2,4] >= -2
 Partition_0_-_Star_edge_{2,_4}: - 2 partition[2,0]
   - circuit_partition_node[4,0] + star_edge[2,4] >= -2
 Partition_0_-_Star_edge_{2,_3}: - 2 partition[3,0]
   - circuit_partition_node[2,0] + star_edge[2,3] >= -2
 Partition_0_-_Star_edge_{2,_3}: - 2 partition[2,0]
   - circuit_partition_node[3,0] + star_edge[2,3] >= -2
 Partition_0_-_Star_edge_{6,_7}: - 2 partition[7,0]
   - circuit_partition_node[6,0] + star_edge[6,7] >= -2
 Partition_0_-_Star_edge_{6,_7}: - 2 partition[6,0]
   - circuit_partition_node[7,0] + star_edge[6,7] >= -2
 Partition_0_-_Star_edge_{6,_8}: - 2 partition[8,0]
   - circuit_partition_node[6,0] + star_edge[6,8] >= -2
 Partition_0_-_Star_edge_{6,_8}: - 2 partition[6,0]
   - circuit_partition_node[8,0] + star_edge[6,8] >= -2
 Partition_0_-_Star_edge_{6,_9}: - 2 partition[9,0]
   - circuit_partition_node[6,0] + star_edge[6,9] >= -2
 Partition_0_-_Star_edge_{6,_9}: - 2 partition[6,0]
   - circuit_partition_node[9,0] + star_edge[6,9] >= -2
 Partition_0_-_Star_edge_{1,_9}: - 2 partition[9,0]
   - circuit_partition_node[1,0] + star_edge[1,9] >= -2
 Partition_0_-_Star_edge_{1,_9}: - 2 partition[1,0]
   - circuit_partition_node[9,0] + star_edge[1,9] >= -2
 Partition_0_-_Star_edge_{1,_8}: - 2 partition[8,0]
   - circuit_partition_node[1,0] + star_edge[1,8] >= -2
 Partition_0_-_Star_edge_{1,_8}: - 2 partition[1,0]
   - circuit_partition_node[8,0] + star_edge[1,8] >= -2
 Partition_0_-_Star_edge_{1,_7}: - 2 partition[7,0]
   - circuit_partition_node[1,0] + star_edge[1,7] >= -2
 Partition_0_-_Star_edge_{1,_7}: - 2 partition[1,0]
   - circuit_partition_node[7,0] + star_edge[1,7] >= -2
 Partition_0_-_Star_edge_{1,_6}: - 2 partition[6,0]
   - circuit_partition_node[1,0] + star_edge[1,6] >= -2
 Partition_0_-_Star_edge_{1,_6}: - 2 partition[1,0]
   - circuit_partition_node[6,0] + star_edge[1,6] >= -2
 Partition_0_-_Star_edge_{1,_5}: - 2 partition[5,0]
   - circuit_partition_node[1,0] + star_edge[1,5] >= -2
 Partition_0_-_Star_edge_{1,_5}: - 2 partition[1,0]
   - circuit_partition_node[5,0] + star_edge[1,5] >= -2
 Partition_0_-_Star_edge_{1,_4}: - 2 partition[4,0]
   - circuit_partition_node[1,0] + star_edge[1,4] >= -2
 Partition_0_-_Star_edge_{1,_4}: - 2 partition[1,0]
   - circuit_partition_node[4,0] + star_edge[1,4] >= -2
 Partition_0_-_Star_edge_{1,_3}: - 2 partition[3,0]
   - circuit_partition_node[1,0] + star_edge[1,3] >= -2
 Partition_0_-_Star_edge_{1,_3}: - 2 partition[1,0]
   - circuit_partition_node[3,0] + star_edge[1,3] >= -2
 Partition_0_-_Star_edge_{1,_2}: - 2 partition[2,0]
   - circuit_partition_node[1,0] + star_edge[1,2] >= -2
 Partition_0_-_Star_edge_{1,_2}: - 2 partition[1,0]
   - circuit_partition_node[2,0] + star_edge[1,2] >= -2
 Partition_0_-_Star_edge_{7,_8}: - 2 partition[8,0]
   - circuit_partition_node[7,0] + star_edge[7,8] >= -2
 Partition_0_-_Star_edge_{7,_8}: - 2 partition[7,0]
   - circuit_partition_node[8,0] + star_edge[7,8] >= -2
 Partition_0_-_Star_edge_{7,_9}: - 2 partition[9,0]
   - circuit_partition_node[7,0] + star_edge[7,9] >= -2
 Partition_0_-_Star_edge_{7,_9}: - 2 partition[7,0]
   - circuit_partition_node[9,0] + star_edge[7,9] >= -2
 Partition_0_-_Star_edge_{0,_9}: - 2 partition[9,0]
   - circuit_partition_node[0,0] + star_edge[0,9] >= -2
 Partition_0_-_Star_edge_{0,_9}: - 2 partition[0,0]
   - circuit_partition_node[9,0] + star_edge[0,9] >= -2
 Partition_0_-_Star_edge_{0,_8}: - 2 partition[8,0]
   - circuit_partition_node[0,0] + star_edge[0,8] >= -2
 Partition_0_-_Star_edge_{0,_8}: - 2 partition[0,0]
   - circuit_partition_node[8,0] + star_edge[0,8] >= -2
 Partition_0_-_Star_edge_{0,_7}: - 2 partition[7,0]
   - circuit_partition_node[0,0] + star_edge[0,7] >= -2
 Partition_0_-_Star_edge_{0,_7}: - 2 partition[0,0]
   - circuit_partition_node[7,0] + star_edge[0,7] >= -2
 Partition_0_-_Star_edge_{0,_6}: - 2 partition[6,0]
   - circuit_partition_node[0,0] + star_edge[0,6] >= -2
 Partition_0_-_Star_edge_{0,_6}: - 2 partition[0,0]
   - circuit_partition_node[6,0] + star_edge[0,6] >= -2
 Partition_0_-_Star_edge_{0,_5}: - 2 partition[5,0]
   - circuit_partition_node[0,0] + star_edge[0,5] >= -2
 Partition_0_-_Star_edge_{0,_5}: - 2 partition[0,0]
   - circuit_partition_node[5,0] + star_edge[0,5] >= -2
 Partition_0_-_Star_edge_{0,_4}: - 2 partition[4,0]
   - circuit_partition_node[0,0] + star_edge[0,4] >= -2
 Partition_0_-_Star_edge_{0,_4}: - 2 partition[0,0]
   - circuit_partition_node[4,0] + star_edge[0,4] >= -2
 Partition_0_-_Star_edge_{0,_3}: - 2 partition[3,0]
   - circuit_partition_node[0,0] + star_edge[0,3] >= -2
 Partition_0_-_Star_edge_{0,_3}: - 2 partition[0,0]
   - circuit_partition_node[3,0] + star_edge[0,3] >= -2
 Partition_0_-_Star_edge_{0,_2}: - 2 partition[2,0]
   - circuit_partition_node[0,0] + star_edge[0,2] >= -2
 Partition_0_-_Star_edge_{0,_2}: - 2 partition[0,0]
   - circuit_partition_node[2,0] + star_edge[0,2] >= -2
 Partition_0_-_Star_edge_{0,_1}: - 2 partition[1,0]
   - circuit_partition_node[0,0] + star_edge[0,1] >= -2
 Partition_0_-_Star_edge_{0,_1}: - 2 partition[0,0]
   - circuit_partition_node[1,0] + star_edge[0,1] >= -2
 Partition_0_-_Star_edge_{8,_9}: - 2 partition[9,0]
   - circuit_partition_node[8,0] + star_edge[8,9] >= -2
 Partition_0_-_Star_edge_{8,_9}: - 2 partition[8,0]
   - circuit_partition_node[9,0] + star_edge[8,9] >= -2
 Partition_1_-_Star_edge_{4,_5}: - 2 partition[5,1]
   - circuit_partition_node[4,1] + star_edge[4,5] >= -2
 Partition_1_-_Star_edge_{4,_5}: - 2 partition[4,1]
   - circuit_partition_node[5,1] + star_edge[4,5] >= -2
 Partition_1_-_Star_edge_{4,_6}: - 2 partition[6,1]
   - circuit_partition_node[4,1] + star_edge[4,6] >= -2
 Partition_1_-_Star_edge_{4,_6}: - 2 partition[4,1]
   - circuit_partition_node[6,1] + star_edge[4,6] >= -2
 Partition_1_-_Star_edge_{4,_7}: - 2 partition[7,1]
   - circuit_partition_node[4,1] + star_edge[4,7] >= -2
 Partition_1_-_Star_edge_{4,_7}: - 2 partition[4,1]
   - circuit_partition_node[7,1] + star_edge[4,7] >= -2
 Partition_1_-_Star_edge_{4,_8}: - 2 partition[8,1]
   - circuit_partition_node[4,1] + star_edge[4,8] >= -2
 Partition_1_-_Star_edge_{4,_8}: - 2 partition[4,1]
   - circuit_partition_node[8,1] + star_edge[4,8] >= -2
 Partition_1_-_Star_edge_{4,_9}: - 2 partition[9,1]
   - circuit_partition_node[4,1] + star_edge[4,9] >= -2
 Partition_1_-_Star_edge_{4,_9}: - 2 partition[4,1]
   - circuit_partition_node[9,1] + star_edge[4,9] >= -2
 Partition_1_-_Star_edge_{3,_9}: - 2 partition[9,1]
   - circuit_partition_node[3,1] + star_edge[3,9] >= -2
 Partition_1_-_Star_edge_{3,_9}: - 2 partition[3,1]
   - circuit_partition_node[9,1] + star_edge[3,9] >= -2
 Partition_1_-_Star_edge_{3,_8}: - 2 partition[8,1]
   - circuit_partition_node[3,1] + star_edge[3,8] >= -2
 Partition_1_-_Star_edge_{3,_8}: - 2 partition[3,1]
   - circuit_partition_node[8,1] + star_edge[3,8] >= -2
 Partition_1_-_Star_edge_{3,_7}: - 2 partition[7,1]
   - circuit_partition_node[3,1] + star_edge[3,7] >= -2
 Partition_1_-_Star_edge_{3,_7}: - 2 partition[3,1]
   - circuit_partition_node[7,1] + star_edge[3,7] >= -2
 Partition_1_-_Star_edge_{3,_6}: - 2 partition[6,1]
   - circuit_partition_node[3,1] + star_edge[3,6] >= -2
 Partition_1_-_Star_edge_{3,_6}: - 2 partition[3,1]
   - circuit_partition_node[6,1] + star_edge[3,6] >= -2
 Partition_1_-_Star_edge_{3,_5}: - 2 partition[5,1]
   - circuit_partition_node[3,1] + star_edge[3,5] >= -2
 Partition_1_-_Star_edge_{3,_5}: - 2 partition[3,1]
   - circuit_partition_node[5,1] + star_edge[3,5] >= -2
 Partition_1_-_Star_edge_{3,_4}: - 2 partition[4,1]
   - circuit_partition_node[3,1] + star_edge[3,4] >= -2
 Partition_1_-_Star_edge_{3,_4}: - 2 partition[3,1]
   - circuit_partition_node[4,1] + star_edge[3,4] >= -2
 Partition_1_-_Star_edge_{5,_6}: - 2 partition[6,1]
   - circuit_partition_node[5,1] + star_edge[5,6] >= -2
 Partition_1_-_Star_edge_{5,_6}: - 2 partition[5,1]
   - circuit_partition_node[6,1] + star_edge[5,6] >= -2
 Partition_1_-_Star_edge_{5,_7}: - 2 partition[7,1]
   - circuit_partition_node[5,1] + star_edge[5,7] >= -2
 Partition_1_-_Star_edge_{5,_7}: - 2 partition[5,1]
   - circuit_partition_node[7,1] + star_edge[5,7] >= -2
 Partition_1_-_Star_edge_{5,_8}: - 2 partition[8,1]
   - circuit_partition_node[5,1] + star_edge[5,8] >= -2
 Partition_1_-_Star_edge_{5,_8}: - 2 partition[5,1]
   - circuit_partition_node[8,1] + star_edge[5,8] >= -2
 Partition_1_-_Star_edge_{5,_9}: - 2 partition[9,1]
   - circuit_partition_node[5,1] + star_edge[5,9] >= -2
 Partition_1_-_Star_edge_{5,_9}: - 2 partition[5,1]
   - circuit_partition_node[9,1] + star_edge[5,9] >= -2
 Partition_1_-_Star_edge_{2,_9}: - 2 partition[9,1]
   - circuit_partition_node[2,1] + star_edge[2,9] >= -2
 Partition_1_-_Star_edge_{2,_9}: - 2 partition[2,1]
   - circuit_partition_node[9,1] + star_edge[2,9] >= -2
 Partition_1_-_Star_edge_{2,_8}: - 2 partition[8,1]
   - circuit_partition_node[2,1] + star_edge[2,8] >= -2
 Partition_1_-_Star_edge_{2,_8}: - 2 partition[2,1]
   - circuit_partition_node[8,1] + star_edge[2,8] >= -2
 Partition_1_-_Star_edge_{2,_7}: - 2 partition[7,1]
   - circuit_partition_node[2,1] + star_edge[2,7] >= -2
 Partition_1_-_Star_edge_{2,_7}: - 2 partition[2,1]
   - circuit_partition_node[7,1] + star_edge[2,7] >= -2
 Partition_1_-_Star_edge_{2,_6}: - 2 partition[6,1]
   - circuit_partition_node[2,1] + star_edge[2,6] >= -2
 Partition_1_-_Star_edge_{2,_6}: - 2 partition[2,1]
   - circuit_partition_node[6,1] + star_edge[2,6] >= -2
 Partition_1_-_Star_edge_{2,_5}: - 2 partition[5,1]
   - circuit_partition_node[2,1] + star_edge[2,5] >= -2
 Partition_1_-_Star_edge_{2,_5}: - 2 partition[2,1]
   - circuit_partition_node[5,1] + star_edge[2,5] >= -2
 Partition_1_-_Star_edge_{2,_4}: - 2 partition[4,1]
   - circuit_partition_node[2,1] + star_edge[2,4] >= -2
 Partition_1_-_Star_edge_{2,_4}: - 2 partition[2,1]
   - circuit_partition_node[4,1] + star_edge[2,4] >= -2
 Partition_1_-_Star_edge_{2,_3}: - 2 partition[3,1]
   - circuit_partition_node[2,1] + star_edge[2,3] >= -2
 Partition_1_-_Star_edge_{2,_3}: - 2 partition[2,1]
   - circuit_partition_node[3,1] + star_edge[2,3] >= -2
 Partition_1_-_Star_edge_{6,_7}: - 2 partition[7,1]
   - circuit_partition_node[6,1] + star_edge[6,7] >= -2
 Partition_1_-_Star_edge_{6,_7}: - 2 partition[6,1]
   - circuit_partition_node[7,1] + star_edge[6,7] >= -2
 Partition_1_-_Star_edge_{6,_8}: - 2 partition[8,1]
   - circuit_partition_node[6,1] + star_edge[6,8] >= -2
 Partition_1_-_Star_edge_{6,_8}: - 2 partition[6,1]
   - circuit_partition_node[8,1] + star_edge[6,8] >= -2
 Partition_1_-_Star_edge_{6,_9}: - 2 partition[9,1]
   - circuit_partition_node[6,1] + star_edge[6,9] >= -2
 Partition_1_-_Star_edge_{6,_9}: - 2 partition[6,1]
   - circuit_partition_node[9,1] + star_edge[6,9] >= -2
 Partition_1_-_Star_edge_{1,_9}: - 2 partition[9,1]
   - circuit_partition_node[1,1] + star_edge[1,9] >= -2
 Partition_1_-_Star_edge_{1,_9}: - 2 partition[1,1]
   - circuit_partition_node[9,1] + star_edge[1,9] >= -2
 Partition_1_-_Star_edge_{1,_8}: - 2 partition[8,1]
   - circuit_partition_node[1,1] + star_edge[1,8] >= -2
 Partition_1_-_Star_edge_{1,_8}: - 2 partition[1,1]
   - circuit_partition_node[8,1] + star_edge[1,8] >= -2
 Partition_1_-_Star_edge_{1,_7}: - 2 partition[7,1]
   - circuit_partition_node[1,1] + star_edge[1,7] >= -2
 Partition_1_-_Star_edge_{1,_7}: - 2 partition[1,1]
   - circuit_partition_node[7,1] + star_edge[1,7] >= -2
 Partition_1_-_Star_edge_{1,_6}: - 2 partition[6,1]
   - circuit_partition_node[1,1] + star_edge[1,6] >= -2
 Partition_1_-_Star_edge_{1,_6}: - 2 partition[1,1]
   - circuit_partition_node[6,1] + star_edge[1,6] >= -2
 Partition_1_-_Star_edge_{1,_5}: - 2 partition[5,1]
   - circuit_partition_node[1,1] + star_edge[1,5] >= -2
 Partition_1_-_Star_edge_{1,_5}: - 2 partition[1,1]
   - circuit_partition_node[5,1] + star_edge[1,5] >= -2
 Partition_1_-_Star_edge_{1,_4}: - 2 partition[4,1]
   - circuit_partition_node[1,1] + star_edge[1,4] >= -2
 Partition_1_-_Star_edge_{1,_4}: - 2 partition[1,1]
   - circuit_partition_node[4,1] + star_edge[1,4] >= -2
 Partition_1_-_Star_edge_{1,_3}: - 2 partition[3,1]
   - circuit_partition_node[1,1] + star_edge[1,3] >= -2
 Partition_1_-_Star_edge_{1,_3}: - 2 partition[1,1]
   - circuit_partition_node[3,1] + star_edge[1,3] >= -2
 Partition_1_-_Star_edge_{1,_2}: - 2 partition[2,1]
   - circuit_partition_node[1,1] + star_edge[1,2] >= -2
 Partition_1_-_Star_edge_{1,_2}: - 2 partition[1,1]
   - circuit_partition_node[2,1] + star_edge[1,2] >= -2
 Partition_1_-_Star_edge_{7,_8}: - 2 partition[8,1]
   - circuit_partition_node[7,1] + star_edge[7,8] >= -2
 Partition_1_-_Star_edge_{7,_8}: - 2 partition[7,1]
   - circuit_partition_node[8,1] + star_edge[7,8] >= -2
 Partition_1_-_Star_edge_{7,_9}: - 2 partition[9,1]
   - circuit_partition_node[7,1] + star_edge[7,9] >= -2
 Partition_1_-_Star_edge_{7,_9}: - 2 partition[7,1]
   - circuit_partition_node[9,1] + star_edge[7,9] >= -2
 Partition_1_-_Star_edge_{0,_9}: - 2 partition[9,1]
   - circuit_partition_node[0,1] + star_edge[0,9] >= -2
 Partition_1_-_Star_edge_{0,_9}: - 2 partition[0,1]
   - circuit_partition_node[9,1] + star_edge[0,9] >= -2
 Partition_1_-_Star_edge_{0,_8}: - 2 partition[8,1]
   - circuit_partition_node[0,1] + star_edge[0,8] >= -2
 Partition_1_-_Star_edge_{0,_8}: - 2 partition[0,1]
   - circuit_partition_node[8,1] + star_edge[0,8] >= -2
 Partition_1_-_Star_edge_{0,_7}: - 2 partition[7,1]
   - circuit_partition_node[0,1] + star_edge[0,7] >= -2
 Partition_1_-_Star_edge_{0,_7}: - 2 partition[0,1]
   - circuit_partition_node[7,1] + star_edge[0,7] >= -2
 Partition_1_-_Star_edge_{0,_6}: - 2 partition[6,1]
   - circuit_partition_node[0,1] + star_edge[0,6] >= -2
 Partition_1_-_Star_edge_{0,_6}: - 2 partition[0,1]
   - circuit_partition_node[6,1] + star_edge[0,6] >= -2
 Partition_1_-_Star_edge_{0,_5}: - 2 partition[5,1]
   - circuit_partition_node[0,1] + star_edge[0,5] >= -2
 Partition_1_-_Star_edge_{0,_5}: - 2 partition[0,1]
   - circuit_partition_node[5,1] + star_edge[0,5] >= -2
 Partition_1_-_Star_edge_{0,_4}: - 2 partition[4,1]
   - circuit_partition_node[0,1] + star_edge[0,4] >= -2
 Partition_1_-_Star_edge_{0,_4}: - 2 partition[0,1]
   - circuit_partition_node[4,1] + star_edge[0,4] >= -2
 Partition_1_-_Star_edge_{0,_3}: - 2 partition[3,1]
   - circuit_partition_node[0,1] + star_edge[0,3] >= -2
 Partition_1_-_Star_edge_{0,_3}: - 2 partition[0,1]
   - circuit_partition_node[3,1] + star_edge[0,3] >= -2
 Partition_1_-_Star_edge_{0,_2}: - 2 partition[2,1]
   - circuit_partition_node[0,1] + star_edge[0,2] >= -2
 Partition_1_-_Star_edge_{0,_2}: - 2 partition[0,1]
   - circuit_partition_node[2,1] + star_edge[0,2] >= -2
 Partition_1_-_Star_edge_{0,_1}: - 2 partition[1,1]
   - circuit_partition_node[0,1] + star_edge[0,1] >= -2
 Partition_1_-_Star_edge_{0,_1}: - 2 partition[0,1]
   - circuit_partition_node[1,1] + star_edge[0,1] >= -2
 Partition_1_-_Star_edge_{8,_9}: - 2 partition[9,1]
   - circuit_partition_node[8,1] + star_edge[8,9] >= -2
 Partition_1_-_Star_edge_{8,_9}: - 2 partition[8,1]
   - circuit_partition_node[9,1] + star_edge[8,9] >= -2
 Partition_2_-_Star_edge_{4,_5}: - 2 partition[5,2]
   - circuit_partition_node[4,2] + star_edge[4,5] >= -2
 Partition_2_-_Star_edge_{4,_5}: - 2 partition[4,2]
   - circuit_partition_node[5,2] + star_edge[4,5] >= -2
 Partition_2_-_Star_edge_{4,_6}: - 2 partition[6,2]
   - circuit_partition_node[4,2] + star_edge[4,6] >= -2
 Partition_2_-_Star_edge_{4,_6}: - 2 partition[4,2]
   - circuit_partition_node[6,2] + star_edge[4,6] >= -2
 Partition_2_-_Star_edge_{4,_7}: - 2 partition[7,2]
   - circuit_partition_node[4,2] + star_edge[4,7] >= -2
 Partition_2_-_Star_edge_{4,_7}: - 2 partition[4,2]
   - circuit_partition_node[7,2] + star_edge[4,7] >= -2
 Partition_2_-_Star_edge_{4,_8}: - 2 partition[8,2]
   - circuit_partition_node[4,2] + star_edge[4,8] >= -2
 Partition_2_-_Star_edge_{4,_8}: - 2 partition[4,2]
   - circuit_partition_node[8,2] + star_edge[4,8] >= -2
 Partition_2_-_Star_edge_{4,_9}: - 2 partition[9,2]
   - circuit_partition_node[4,2] + star_edge[4,9] >= -2
 Partition_2_-_Star_edge_{4,_9}: - 2 partition[4,2]
   - circuit_partition_node[9,2] + star_edge[4,9] >= -2
 Partition_2_-_Star_edge_{3,_9}: - 2 partition[9,2]
   - circuit_partition_node[3,2] + star_edge[3,9] >= -2
 Partition_2_-_Star_edge_{3,_9}: - 2 partition[3,2]
   - circuit_partition_node[9,2] + star_edge[3,9] >= -2
 Partition_2_-_Star_edge_{3,_8}: - 2 partition[8,2]
   - circuit_partition_node[3,2] + star_edge[3,8] >= -2
 Partition_2_-_Star_edge_{3,_8}: - 2 partition[3,2]
   - circuit_partition_node[8,2] + star_edge[3,8] >= -2
 Partition_2_-_Star_edge_{3,_7}: - 2 partition[7,2]
   - circuit_partition_node[3,2] + star_edge[3,7] >= -2
 Partition_2_-_Star_edge_{3,_7}: - 2 partition[3,2]
   - circuit_partition_node[7,2] + star_edge[3,7] >= -2
 Partition_2_-_Star_edge_{3,_6}: - 2 partition[6,2]
   - circuit_partition_node[3,2] + star_edge[3,6] >= -2
 Partition_2_-_Star_edge_{3,_6}: - 2 partition[3,2]
   - circuit_partition_node[6,2] + star_edge[3,6] >= -2
 Partition_2_-_Star_edge_{3,_5}: - 2 partition[5,2]
   - circuit_partition_node[3,2] + star_edge[3,5] >= -2
 Partition_2_-_Star_edge_{3,_5}: - 2 partition[3,2]
   - circuit_partition_node[5,2] + star_edge[3,5] >= -2
 Partition_2_-_Star_edge_{3,_4}: - 2 partition[4,2]
   - circuit_partition_node[3,2] + star_edge[3,4] >= -2
 Partition_2_-_Star_edge_{3,_4}: - 2 partition[3,2]
   - circuit_partition_node[4,2] + star_edge[3,4] >= -2
 Partition_2_-_Star_edge_{5,_6}: - 2 partition[6,2]
   - circuit_partition_node[5,2] + star_edge[5,6] >= -2
 Partition_2_-_Star_edge_{5,_6}: - 2 partition[5,2]
   - circuit_partition_node[6,2] + star_edge[5,6] >= -2
 Partition_2_-_Star_edge_{5,_7}: - 2 partition[7,2]
   - circuit_partition_node[5,2] + star_edge[5,7] >= -2
 Partition_2_-_Star_edge_{5,_7}: - 2 partition[5,2]
   - circuit_partition_node[7,2] + star_edge[5,7] >= -2
 Partition_2_-_Star_edge_{5,_8}: - 2 partition[8,2]
   - circuit_partition_node[5,2] + star_edge[5,8] >= -2
 Partition_2_-_Star_edge_{5,_8}: - 2 partition[5,2]
   - circuit_partition_node[8,2] + star_edge[5,8] >= -2
 Partition_2_-_Star_edge_{5,_9}: - 2 partition[9,2]
   - circuit_partition_node[5,2] + star_edge[5,9] >= -2
 Partition_2_-_Star_edge_{5,_9}: - 2 partition[5,2]
   - circuit_partition_node[9,2] + star_edge[5,9] >= -2
 Partition_2_-_Star_edge_{2,_9}: - 2 partition[9,2]
   - circuit_partition_node[2,2] + star_edge[2,9] >= -2
 Partition_2_-_Star_edge_{2,_9}: - 2 partition[2,2]
   - circuit_partition_node[9,2] + star_edge[2,9] >= -2
 Partition_2_-_Star_edge_{2,_8}: - 2 partition[8,2]
   - circuit_partition_node[2,2] + star_edge[2,8] >= -2
 Partition_2_-_Star_edge_{2,_8}: - 2 partition[2,2]
   - circuit_partition_node[8,2] + star_edge[2,8] >= -2
 Partition_2_-_Star_edge_{2,_7}: - 2 partition[7,2]
   - circuit_partition_node[2,2] + star_edge[2,7] >= -2
 Partition_2_-_Star_edge_{2,_7}: - 2 partition[2,2]
   - circuit_partition_node[7,2] + star_edge[2,7] >= -2
 Partition_2_-_Star_edge_{2,_6}: - 2 partition[6,2]
   - circuit_partition_node[2,2] + star_edge[2,6] >= -2
 Partition_2_-_Star_edge_{2,_6}: - 2 partition[2,2]
   - circuit_partition_node[6,2] + star_edge[2,6] >= -2
 Partition_2_-_Star_edge_{2,_5}: - 2 partition[5,2]
   - circuit_partition_node[2,2] + star_edge[2,5] >= -2
 Partition_2_-_Star_edge_{2,_5}: - 2 partition[2,2]
   - circuit_partition_node[5,2] + star_edge[2,5] >= -2
 Partition_2_-_Star_edge_{2,_4}: - 2 partition[4,2]
   - circuit_partition_node[2,2] + star_edge[2,4] >= -2
 Partition_2_-_Star_edge_{2,_4}: - 2 partition[2,2]
   - circuit_partition_node[4,2] + star_edge[2,4] >= -2
 Partition_2_-_Star_edge_{2,_3}: - 2 partition[3,2]
   - circuit_partition_node[2,2] + star_edge[2,3] >= -2
 Partition_2_-_Star_edge_{2,_3}: - 2 partition[2,2]
   - circuit_partition_node[3,2] + star_edge[2,3] >= -2
 Partition_2_-_Star_edge_{6,_7}: - 2 partition[7,2]
   - circuit_partition_node[6,2] + star_edge[6,7] >= -2
 Partition_2_-_Star_edge_{6,_7}: - 2 partition[6,2]
   - circuit_partition_node[7,2] + star_edge[6,7] >= -2
 Partition_2_-_Star_edge_{6,_8}: - 2 partition[8,2]
   - circuit_partition_node[6,2] + star_edge[6,8] >= -2
 Partition_2_-_Star_edge_{6,_8}: - 2 partition[6,2]
   - circuit_partition_node[8,2] + star_edge[6,8] >= -2
 Partition_2_-_Star_edge_{6,_9}: - 2 partition[9,2]
   - circuit_partition_node[6,2] + star_edge[6,9] >= -2
 Partition_2_-_Star_edge_{6,_9}: - 2 partition[6,2]
   - circuit_partition_node[9,2] + star_edge[6,9] >= -2
 Partition_2_-_Star_edge_{1,_9}: - 2 partition[9,2]
   - circuit_partition_node[1,2] + star_edge[1,9] >= -2
 Partition_2_-_Star_edge_{1,_9}: - 2 partition[1,2]
   - circuit_partition_node[9,2] + star_edge[1,9] >= -2
 Partition_2_-_Star_edge_{1,_8}: - 2 partition[8,2]
   - circuit_partition_node[1,2] + star_edge[1,8] >= -2
 Partition_2_-_Star_edge_{1,_8}: - 2 partition[1,2]
   - circuit_partition_node[8,2] + star_edge[1,8] >= -2
 Partition_2_-_Star_edge_{1,_7}: - 2 partition[7,2]
   - circuit_partition_node[1,2] + star_edge[1,7] >= -2
 Partition_2_-_Star_edge_{1,_7}: - 2 partition[1,2]
   - circuit_partition_node[7,2] + star_edge[1,7] >= -2
 Partition_2_-_Star_edge_{1,_6}: - 2 partition[6,2]
   - circuit_partition_node[1,2] + star_edge[1,6] >= -2
 Partition_2_-_Star_edge_{1,_6}: - 2 partition[1,2]
   - circuit_partition_node[6,2] + star_edge[1,6] >= -2
 Partition_2_-_Star_edge_{1,_5}: - 2 partition[5,2]
   - circuit_partition_node[1,2] + star_edge[1,5] >= -2
 Partition_2_-_Star_edge_{1,_5}: - 2 partition[1,2]
   - circuit_partition_node[5,2] + star_edge[1,5] >= -2
 Partition_2_-_Star_edge_{1,_4}: - 2 partition[4,2]
   - circuit_partition_node[1,2] + star_edge[1,4] >= -2
 Partition_2_-_Star_edge_{1,_4}: - 2 partition[1,2]
   - circuit_partition_node[4,2] + star_edge[1,4] >= -2
 Partition_2_-_Star_edge_{1,_3}: - 2 partition[3,2]
   - circuit_partition_node[1,2] + star_edge[1,3] >= -2
 Partition_2_-_Star_edge_{1,_3}: - 2 partition[1,2]
   - circuit_partition_node[3,2] + star_edge[1,3] >= -2
 Partition_2_-_Star_edge_{1,_2}: - 2 partition[2,2]
   - circuit_partition_node[1,2] + star_edge[1,2] >= -2
 Partition_2_-_Star_edge_{1,_2}: - 2 partition[1,2]
   - circuit_partition_node[2,2] + star_edge[1,2] >= -2
 Partition_2_-_Star_edge_{7,_8}: - 2 partition[8,2]
   - circuit_partition_node[7,2] + star_edge[7,8] >= -2
 Partition_2_-_Star_edge_{7,_8}: - 2 partition[7,2]
   - circuit_partition_node[8,2] + star_edge[7,8] >= -2
 Partition_2_-_Star_edge_{7,_9}: - 2 partition[9,2]
   - circuit_partition_node[7,2] + star_edge[7,9] >= -2
 Partition_2_-_Star_edge_{7,_9}: - 2 partition[7,2]
   - circuit_partition_node[9,2] + star_edge[7,9] >= -2
 Partition_2_-_Star_edge_{0,_9}: - 2 partition[9,2]
   - circuit_partition_node[0,2] + star_edge[0,9] >= -2
 Partition_2_-_Star_edge_{0,_9}: - 2 partition[0,2]
   - circuit_partition_node[9,2] + star_edge[0,9] >= -2
 Partition_2_-_Star_edge_{0,_8}: - 2 partition[8,2]
   - circuit_partition_node[0,2] + star_edge[0,8] >= -2
 Partition_2_-_Star_edge_{0,_8}: - 2 partition[0,2]
   - circuit_partition_node[8,2] + star_edge[0,8] >= -2
 Partition_2_-_Star_edge_{0,_7}: - 2 partition[7,2]
   - circuit_partition_node[0,2] + star_edge[0,7] >= -2
 Partition_2_-_Star_edge_{0,_7}: - 2 partition[0,2]
   - circuit_partition_node[7,2] + star_edge[0,7] >= -2
 Partition_2_-_Star_edge_{0,_6}: - 2 partition[6,2]
   - circuit_partition_node[0,2] + star_edge[0,6] >= -2
 Partition_2_-_Star_edge_{0,_6}: - 2 partition[0,2]
   - circuit_partition_node[6,2] + star_edge[0,6] >= -2
 Partition_2_-_Star_edge_{0,_5}: - 2 partition[5,2]
   - circuit_partition_node[0,2] + star_edge[0,5] >= -2
 Partition_2_-_Star_edge_{0,_5}: - 2 partition[0,2]
   - circuit_partition_node[5,2] + star_edge[0,5] >= -2
 Partition_2_-_Star_edge_{0,_4}: - 2 partition[4,2]
   - circuit_partition_node[0,2] + star_edge[0,4] >= -2
 Partition_2_-_Star_edge_{0,_4}: - 2 partition[0,2]
   - circuit_partition_node[4,2] + star_edge[0,4] >= -2
 Partition_2_-_Star_edge_{0,_3}: - 2 partition[3,2]
   - circuit_partition_node[0,2] + star_edge[0,3] >= -2
 Partition_2_-_Star_edge_{0,_3}: - 2 partition[0,2]
   - circuit_partition_node[3,2] + star_edge[0,3] >= -2
 Partition_2_-_Star_edge_{0,_2}: - 2 partition[2,2]
   - circuit_partition_node[0,2] + star_edge[0,2] >= -2
 Partition_2_-_Star_edge_{0,_2}: - 2 partition[0,2]
   - circuit_partition_node[2,2] + star_edge[0,2] >= -2
 Partition_2_-_Star_edge_{0,_1}: - 2 partition[1,2]
   - circuit_partition_node[0,2] + star_edge[0,1] >= -2
 Partition_2_-_Star_edge_{0,_1}: - 2 partition[0,2]
   - circuit_partition_node[1,2] + star_edge[0,1] >= -2
 Partition_2_-_Star_edge_{8,_9}: - 2 partition[9,2]
   - circuit_partition_node[8,2] + star_edge[8,9] >= -2
 Partition_2_-_Star_edge_{8,_9}: - 2 partition[8,2]
   - circuit_partition_node[9,2] + star_edge[8,9] >= -2
 Partition_3_-_Star_edge_{4,_5}: - 2 partition[5,3]
   - circuit_partition_node[4,3] + star_edge[4,5] >= -2
 Partition_3_-_Star_edge_{4,_5}: - 2 partition[4,3]
   - circuit_partition_node[5,3] + star_edge[4,5] >= -2
 Partition_3_-_Star_edge_{4,_6}: - 2 partition[6,3]
   - circuit_partition_node[4,3] + star_edge[4,6] >= -2
 Partition_3_-_Star_edge_{4,_6}: - 2 partition[4,3]
   - circuit_partition_node[6,3] + star_edge[4,6] >= -2
 Partition_3_-_Star_edge_{4,_7}: - 2 partition[7,3]
   - circuit_partition_node[4,3] + star_edge[4,7] >= -2
 Partition_3_-_Star_edge_{4,_7}: - 2 partition[4,3]
   - circuit_partition_node[7,3] + star_edge[4,7] >= -2
 Partition_3_-_Star_edge_{4,_8}: - 2 partition[8,3]
   - circuit_partition_node[4,3] + star_edge[4,8] >= -2
 Partition_3_-_Star_edge_{4,_8}: - 2 partition[4,3]
   - circuit_partition_node[8,3] + star_edge[4,8] >= -2
 Partition_3_-_Star_edge_{4,_9}: - 2 partition[9,3]
   - circuit_partition_node[4,3] + star_edge[4,9] >= -2
 Partition_3_-_Star_edge_{4,_9}: - 2 partition[4,3]
   - circuit_partition_node[9,3] + star_edge[4,9] >= -2
 Partition_3_-_Star_edge_{3,_9}: - 2 partition[9,3]
   - circuit_partition_node[3,3] + star_edge[3,9] >= -2
 Partition_3_-_Star_edge_{3,_9}: - 2 partition[3,3]
   - circuit_partition_node[9,3] + star_edge[3,9] >= -2
 Partition_3_-_Star_edge_{3,_8}: - 2 partition[8,3]
   - circuit_partition_node[3,3] + star_edge[3,8] >= -2
 Partition_3_-_Star_edge_{3,_8}: - 2 partition[3,3]
   - circuit_partition_node[8,3] + star_edge[3,8] >= -2
 Partition_3_-_Star_edge_{3,_7}: - 2 partition[7,3]
   - circuit_partition_node[3,3] + star_edge[3,7] >= -2
 Partition_3_-_Star_edge_{3,_7}: - 2 partition[3,3]
   - circuit_partition_node[7,3] + star_edge[3,7] >= -2
 Partition_3_-_Star_edge_{3,_6}: - 2 partition[6,3]
   - circuit_partition_node[3,3] + star_edge[3,6] >= -2
 Partition_3_-_Star_edge_{3,_6}: - 2 partition[3,3]
   - circuit_partition_node[6,3] + star_edge[3,6] >= -2
 Partition_3_-_Star_edge_{3,_5}: - 2 partition[5,3]
   - circuit_partition_node[3,3] + star_edge[3,5] >= -2
 Partition_3_-_Star_edge_{3,_5}: - 2 partition[3,3]
   - circuit_partition_node[5,3] + star_edge[3,5] >= -2
 Partition_3_-_Star_edge_{3,_4}: - 2 partition[4,3]
   - circuit_partition_node[3,3] + star_edge[3,4] >= -2
 Partition_3_-_Star_edge_{3,_4}: - 2 partition[3,3]
   - circuit_partition_node[4,3] + star_edge[3,4] >= -2
 Partition_3_-_Star_edge_{5,_6}: - 2 partition[6,3]
   - circuit_partition_node[5,3] + star_edge[5,6] >= -2
 Partition_3_-_Star_edge_{5,_6}: - 2 partition[5,3]
   - circuit_partition_node[6,3] + star_edge[5,6] >= -2
 Partition_3_-_Star_edge_{5,_7}: - 2 partition[7,3]
   - circuit_partition_node[5,3] + star_edge[5,7] >= -2
 Partition_3_-_Star_edge_{5,_7}: - 2 partition[5,3]
   - circuit_partition_node[7,3] + star_edge[5,7] >= -2
 Partition_3_-_Star_edge_{5,_8}: - 2 partition[8,3]
   - circuit_partition_node[5,3] + star_edge[5,8] >= -2
 Partition_3_-_Star_edge_{5,_8}: - 2 partition[5,3]
   - circuit_partition_node[8,3] + star_edge[5,8] >= -2
 Partition_3_-_Star_edge_{5,_9}: - 2 partition[9,3]
   - circuit_partition_node[5,3] + star_edge[5,9] >= -2
 Partition_3_-_Star_edge_{5,_9}: - 2 partition[5,3]
   - circuit_partition_node[9,3] + star_edge[5,9] >= -2
 Partition_3_-_Star_edge_{2,_9}: - 2 partition[9,3]
   - circuit_partition_node[2,3] + star_edge[2,9] >= -2
 Partition_3_-_Star_edge_{2,_9}: - 2 partition[2,3]
   - circuit_partition_node[9,3] + star_edge[2,9] >= -2
 Partition_3_-_Star_edge_{2,_8}: - 2 partition[8,3]
   - circuit_partition_node[2,3] + star_edge[2,8] >= -2
 Partition_3_-_Star_edge_{2,_8}: - 2 partition[2,3]
   - circuit_partition_node[8,3] + star_edge[2,8] >= -2
 Partition_3_-_Star_edge_{2,_7}: - 2 partition[7,3]
   - circuit_partition_node[2,3] + star_edge[2,7] >= -2
 Partition_3_-_Star_edge_{2,_7}: - 2 partition[2,3]
   - circuit_partition_node[7,3] + star_edge[2,7] >= -2
 Partition_3_-_Star_edge_{2,_6}: - 2 partition[6,3]
   - circuit_partition_node[2,3] + star_edge[2,6] >= -2
 Partition_3_-_Star_edge_{2,_6}: - 2 partition[2,3]
   - circuit_partition_node[6,3] + star_edge[2,6] >= -2
 Partition_3_-_Star_edge_{2,_5}: - 2 partition[5,3]
   - circuit_partition_node[2,3] + star_edge[2,5] >= -2
 Partition_3_-_Star_edge_{2,_5}: - 2 partition[2,3]
   - circuit_partition_node[5,3] + star_edge[2,5] >= -2
 Partition_3_-_Star_edge_{2,_4}: - 2 partition[4,3]
   - circuit_partition_node[2,3] + star_edge[2,4] >= -2
 Partition_3_-_Star_edge_{2,_4}: - 2 partition[2,3]
   - circuit_partition_node[4,3] + star_edge[2,4] >= -2
 Partition_3_-_Star_edge_{2,_3}: - 2 partition[3,3]
   - circuit_partition_node[2,3] + star_edge[2,3] >= -2
 Partition_3_-_Star_edge_{2,_3}: - 2 partition[2,3]
   - circuit_partition_node[3,3] + star_edge[2,3] >= -2
 Partition_3_-_Star_edge_{6,_7}: - 2 partition[7,3]
   - circuit_partition_node[6,3] + star_edge[6,7] >= -2
 Partition_3_-_Star_edge_{6,_7}: - 2 partition[6,3]
   - circuit_partition_node[7,3] + star_edge[6,7] >= -2
 Partition_3_-_Star_edge_{6,_8}: - 2 partition[8,3]
   - circuit_partition_node[6,3] + star_edge[6,8] >= -2
 Partition_3_-_Star_edge_{6,_8}: - 2 partition[6,3]
   - circuit_partition_node[8,3] + star_edge[6,8] >= -2
 Partition_3_-_Star_edge_{6,_9}: - 2 partition[9,3]
   - circuit_partition_node[6,3] + star_edge[6,9] >= -2
 Partition_3_-_Star_edge_{6,_9}: - 2 partition[6,3]
   - circuit_partition_node[9,3] + star_edge[6,9] >= -2
 Partition_3_-_Star_edge_{1,_9}: - 2 partition[9,3]
   - circuit_partition_node[1,3] + star_edge[1,9] >= -2
 Partition_3_-_Star_edge_{1,_9}: - 2 partition[1,3]
   - circuit_partition_node[9,3] + star_edge[1,9] >= -2
 Partition_3_-_Star_edge_{1,_8}: - 2 partition[8,3]
   - circuit_partition_node[1,3] + star_edge[1,8] >= -2
 Partition_3_-_Star_edge_{1,_8}: - 2 partition[1,3]
   - circuit_partition_node[8,3] + star_edge[1,8] >= -2
 Partition_3_-_Star_edge_{1,_7}: - 2 partition[7,3]
   - circuit_partition_node[1,3] + star_edge[1,7] >= -2
 Partition_3_-_Star_edge_{1,_7}: - 2 partition[1,3]
   - circuit_partition_node[7,3] + star_edge[1,7] >= -2
 Partition_3_-_Star_edge_{1,_6}: - 2 partition[6,3]
   - circuit_partition_node[1,3] + star_edge[1,6] >= -2
 Partition_3_-_Star_edge_{1,_6}: - 2 partition[1,3]
   - circuit_partition_node[6,3] + star_edge[1,6] >= -2
 Partition_3_-_Star_edge_{1,_5}: - 2 partition[5,3]
   - circuit_partition_node[1,3] + star_edge[1,5] >= -2
 Partition_3_-_Star_edge_{1,_5}: - 2 partition[1,3]
   - circuit_partition_node[5,3] + star_edge[1,5] >= -2
 Partition_3_-_Star_edge_{1,_4}: - 2 partition[4,3]
   - circuit_partition_node[1,3] + star_edge[1,4] >= -2
 Partition_3_-_Star_edge_{1,_4}: - 2 partition[1,3]
   - circuit_partition_node[4,3] + star_edge[1,4] >= -2
 Partition_3_-_Star_edge_{1,_3}: - 2 partition[3,3]
   - circuit_partition_node[1,3] + star_edge[1,3] >= -2
 Partition_3_-_Star_edge_{1,_3}: - 2 partition[1,3]
   - circuit_partition_node[3,3] + star_edge[1,3] >= -2
 Partition_3_-_Star_edge_{1,_2}: - 2 partition[2,3]
   - circuit_partition_node[1,3] + star_edge[1,2] >= -2
 Partition_3_-_Star_edge_{1,_2}: - 2 partition[1,3]
   - circuit_partition_node[2,3] + star_edge[1,2] >= -2
 Partition_3_-_Star_edge_{7,_8}: - 2 partition[8,3]
   - circuit_partition_node[7,3] + star_edge[7,8] >= -2
 Partition_3_-_Star_edge_{7,_8}: - 2 partition[7,3]
   - circuit_partition_node[8,3] + star_edge[7,8] >= -2
 Partition_3_-_Star_edge_{7,_9}: - 2 partition[9,3]
   - circuit_partition_node[7,3] + star_edge[7,9] >= -2
 Partition_3_-_Star_edge_{7,_9}: - 2 partition[7,3]
   - circuit_partition_node[9,3] + star_edge[7,9] >= -2
 Partition_3_-_Star_edge_{0,_9}: - 2 partition[9,3]
   - circuit_partition_node[0,3] + star_edge[0,9] >= -2
 Partition_3_-_Star_edge_{0,_9}: - 2 partition[0,3]
   - circuit_partition_node[9,3] + star_edge[0,9] >= -2
 Partition_3_-_Star_edge_{0,_8}: - 2 partition[8,3]
   - circuit_partition_node[0,3] + star_edge[0,8] >= -2
 Partition_3_-_Star_edge_{0,_8}: - 2 partition[0,3]
   - circuit_partition_node[8,3] + star_edge[0,8] >= -2
 Partition_3_-_Star_edge_{0,_7}: - 2 partition[7,3]
   - circuit_partition_node[0,3] + star_edge[0,7] >= -2
 Partition_3_-_Star_edge_{0,_7}: - 2 partition[0,3]
   - circuit_partition_node[7,3] + star_edge[0,7] >= -2
 Partition_3_-_Star_edge_{0,_6}: - 2 partition[6,3]
   - circuit_partition_node[0,3] + star_edge[0,6] >= -2
 Partition_3_-_Star_edge_{0,_6}: - 2 partition[0,3]
   - circuit_partition_node[6,3] + star_edge[0,6] >= -2
 Partition_3_-_Star_edge_{0,_5}: - 2 partition[5,3]
   - circuit_partition_node[0,3] + star_edge[0,5] >= -2
 Partition_3_-_Star_edge_{0,_5}: - 2 partition[0,3]
   - circuit_partition_node[5,3] + star_edge[0,5] >= -2
 Partition_3_-_Star_edge_{0,_4}: - 2 partition[4,3]
   - circuit_partition_node[0,3] + star_edge[0,4] >= -2
 Partition_3_-_Star_edge_{0,_4}: - 2 partition[0,3]
   - circuit_partition_node[4,3] + star_edge[0,4] >= -2
 Partition_3_-_Star_edge_{0,_3}: - 2 partition[3,3]
   - circuit_partition_node[0,3] + star_edge[0,3] >= -2
 Partition_3_-_Star_edge_{0,_3}: - 2 partition[0,3]
   - circuit_partition_node[3,3] + star_edge[0,3] >= -2
 Partition_3_-_Star_edge_{0,_2}: - 2 partition[2,3]
   - circuit_partition_node[0,3] + star_edge[0,2] >= -2
 Partition_3_-_Star_edge_{0,_2}: - 2 partition[0,3]
   - circuit_partition_node[2,3] + star_edge[0,2] >= -2
 Partition_3_-_Star_edge_{0,_1}: - 2 partition[1,3]
   - circuit_partition_node[0,3] + star_edge[0,1] >= -2
 Partition_3_-_Star_edge_{0,_1}: - 2 partition[0,3]
   - circuit_partition_node[1,3] + star_edge[0,1] >= -2
 Partition_3_-_Star_edge_{8,_9}: - 2 partition[9,3]
   - circuit_partition_node[8,3] + star_edge[8,9] >= -2
 Partition_3_-_Star_edge_{8,_9}: - 2 partition[8,3]
   - circuit_partition_node[9,3] + star_edge[8,9] >= -2
 Partition_4_-_Star_edge_{4,_5}: - 2 partition[5,4]
   - circuit_partition_node[4,4] + star_edge[4,5] >= -2
 Partition_4_-_Star_edge_{4,_5}: - 2 partition[4,4]
   - circuit_partition_node[5,4] + star_edge[4,5] >= -2
 Partition_4_-_Star_edge_{4,_6}: - 2 partition[6,4]
   - circuit_partition_node[4,4] + star_edge[4,6] >= -2
 Partition_4_-_Star_edge_{4,_6}: - 2 partition[4,4]
   - circuit_partition_node[6,4] + star_edge[4,6] >= -2
 Partition_4_-_Star_edge_{4,_7}: - 2 partition[7,4]
   - circuit_partition_node[4,4] + star_edge[4,7] >= -2
 Partition_4_-_Star_edge_{4,_7}: - 2 partition[4,4]
   - circuit_partition_node[7,4] + star_edge[4,7] >= -2
 Partition_4_-_Star_edge_{4,_8}: - 2 partition[8,4]
   - circuit_partition_node[4,4] + star_edge[4,8] >= -2
 Partition_4_-_Star_edge_{4,_8}: - 2 partition[4,4]
   - circuit_partition_node[8,4] + star_edge[4,8] >= -2
 Partition_4_-_Star_edge_{4,_9}: - 2 partition[9,4]
   - circuit_partition_node[4,4] + star_edge[4,9] >= -2
 Partition_4_-_Star_edge_{4,_9}: - 2 partition[4,4]
   - circuit_partition_node[9,4] + star_edge[4,9] >= -2
 Partition_4_-_Star_edge_{3,_9}: - 2 partition[9,4]
   - circuit_partition_node[3,4] + star_edge[3,9] >= -2
 Partition_4_-_Star_edge_{3,_9}: - 2 partition[3,4]
   - circuit_partition_node[9,4] + star_edge[3,9] >= -2
 Partition_4_-_Star_edge_{3,_8}: - 2 partition[8,4]
   - circuit_partition_node[3,4] + star_edge[3,8] >= -2
 Partition_4_-_Star_edge_{3,_8}: - 2 partition[3,4]
   - circuit_partition_node[8,4] + star_edge[3,8] >= -2
 Partition_4_-_Star_edge_{3,_7}: - 2 partition[7,4]
   - circuit_partition_node[3,4] + star_edge[3,7] >= -2
 Partition_4_-_Star_edge_{3,_7}: - 2 partition[3,4]
   - circuit_partition_node[7,4] + star_edge[3,7] >= -2
 Partition_4_-_Star_edge_{3,_6}: - 2 partition[6,4]
   - circuit_partition_node[3,4] + star_edge[3,6] >= -2
 Partition_4_-_Star_edge_{3,_6}: - 2 partition[3,4]
   - circuit_partition_node[6,4] + star_edge[3,6] >= -2
 Partition_4_-_Star_edge_{3,_5}: - 2 partition[5,4]
   - circuit_partition_node[3,4] + star_edge[3,5] >= -2
 Partition_4_-_Star_edge_{3,_5}: - 2 partition[3,4]
   - circuit_partition_node[5,4] + star_edge[3,5] >= -2
 Partition_4_-_Star_edge_{3,_4}: - 2 partition[4,4]
   - circuit_partition_node[3,4] + star_edge[3,4] >= -2
 Partition_4_-_Star_edge_{3,_4}: - 2 partition[3,4]
   - circuit_partition_node[4,4] + star_edge[3,4] >= -2
 Partition_4_-_Star_edge_{5,_6}: - 2 partition[6,4]
   - circuit_partition_node[5,4] + star_edge[5,6] >= -2
 Partition_4_-_Star_edge_{5,_6}: - 2 partition[5,4]
   - circuit_partition_node[6,4] + star_edge[5,6] >= -2
 Partition_4_-_Star_edge_{5,_7}: - 2 partition[7,4]
   - circuit_partition_node[5,4] + star_edge[5,7] >= -2
 Partition_4_-_Star_edge_{5,_7}: - 2 partition[5,4]
   - circuit_partition_node[7,4] + star_edge[5,7] >= -2
 Partition_4_-_Star_edge_{5,_8}: - 2 partition[8,4]
   - circuit_partition_node[5,4] + star_edge[5,8] >= -2
 Partition_4_-_Star_edge_{5,_8}: - 2 partition[5,4]
   - circuit_partition_node[8,4] + star_edge[5,8] >= -2
 Partition_4_-_Star_edge_{5,_9}: - 2 partition[9,4]
   - circuit_partition_node[5,4] + star_edge[5,9] >= -2
 Partition_4_-_Star_edge_{5,_9}: - 2 partition[5,4]
   - circuit_partition_node[9,4] + star_edge[5,9] >= -2
 Partition_4_-_Star_edge_{2,_9}: - 2 partition[9,4]
   - circuit_partition_node[2,4] + star_edge[2,9] >= -2
 Partition_4_-_Star_edge_{2,_9}: - 2 partition[2,4]
   - circuit_partition_node[9,4] + star_edge[2,9] >= -2
 Partition_4_-_Star_edge_{2,_8}: - 2 partition[8,4]
   - circuit_partition_node[2,4] + star_edge[2,8] >= -2
 Partition_4_-_Star_edge_{2,_8}: - 2 partition[2,4]
   - circuit_partition_node[8,4] + star_edge[2,8] >= -2
 Partition_4_-_Star_edge_{2,_7}: - 2 partition[7,4]
   - circuit_partition_node[2,4] + star_edge[2,7] >= -2
 Partition_4_-_Star_edge_{2,_7}: - 2 partition[2,4]
   - circuit_partition_node[7,4] + star_edge[2,7] >= -2
 Partition_4_-_Star_edge_{2,_6}: - 2 partition[6,4]
   - circuit_partition_node[2,4] + star_edge[2,6] >= -2
 Partition_4_-_Star_edge_{2,_6}: - 2 partition[2,4]
   - circuit_partition_node[6,4] + star_edge[2,6] >= -2
 Partition_4_-_Star_edge_{2,_5}: - 2 partition[5,4]
   - circuit_partition_node[2,4] + star_edge[2,5] >= -2
 Partition_4_-_Star_edge_{2,_5}: - 2 partition[2,4]
   - circuit_partition_node[5,4] + star_edge[2,5] >= -2
 Partition_4_-_Star_edge_{2,_4}: - 2 partition[4,4]
   - circuit_partition_node[2,4] + star_edge[2,4] >= -2
 Partition_4_-_Star_edge_{2,_4}: - 2 partition[2,4]
   - circuit_partition_node[4,4] + star_edge[2,4] >= -2
 Partition_4_-_Star_edge_{2,_3}: - 2 partition[3,4]
   - circuit_partition_node[2,4] + star_edge[2,3] >= -2
 Partition_4_-_Star_edge_{2,_3}: - 2 partition[2,4]
   - circuit_partition_node[3,4] + star_edge[2,3] >= -2
 Partition_4_-_Star_edge_{6,_7}: - 2 partition[7,4]
   - circuit_partition_node[6,4] + star_edge[6,7] >= -2
 Partition_4_-_Star_edge_{6,_7}: - 2 partition[6,4]
   - circuit_partition_node[7,4] + star_edge[6,7] >= -2
 Partition_4_-_Star_edge_{6,_8}: - 2 partition[8,4]
   - circuit_partition_node[6,4] + star_edge[6,8] >= -2
 Partition_4_-_Star_edge_{6,_8}: - 2 partition[6,4]
   - circuit_partition_node[8,4] + star_edge[6,8] >= -2
 Partition_4_-_Star_edge_{6,_9}: - 2 partition[9,4]
   - circuit_partition_node[6,4] + star_edge[6,9] >= -2
 Partition_4_-_Star_edge_{6,_9}: - 2 partition[6,4]
   - circuit_partition_node[9,4] + star_edge[6,9] >= -2
 Partition_4_-_Star_edge_{1,_9}: - 2 partition[9,4]
   - circuit_partition_node[1,4] + star_edge[1,9] >= -2
 Partition_4_-_Star_edge_{1,_9}: - 2 partition[1,4]
   - circuit_partition_node[9,4] + star_edge[1,9] >= -2
 Partition_4_-_Star_edge_{1,_8}: - 2 partition[8,4]
   - circuit_partition_node[1,4] + star_edge[1,8] >= -2
 Partition_4_-_Star_edge_{1,_8}: - 2 partition[1,4]
   - circuit_partition_node[8,4] + star_edge[1,8] >= -2
 Partition_4_-_Star_edge_{1,_7}: - 2 partition[7,4]
   - circuit_partition_node[1,4] + star_edge[1,7] >= -2
 Partition_4_-_Star_edge_{1,_7}: - 2 partition[1,4]
   - circuit_partition_node[7,4] + star_edge[1,7] >= -2
 Partition_4_-_Star_edge_{1,_6}: - 2 partition[6,4]
   - circuit_partition_node[1,4] + star_edge[1,6] >= -2
 Partition_4_-_Star_edge_{1,_6}: - 2 partition[1,4]
   - circuit_partition_node[6,4] + star_edge[1,6] >= -2
 Partition_4_-_Star_edge_{1,_5}: - 2 partition[5,4]
   - circuit_partition_node[1,4] + star_edge[1,5] >= -2
 Partition_4_-_Star_edge_{1,_5}: - 2 partition[1,4]
   - circuit_partition_node[5,4] + star_edge[1,5] >= -2
 Partition_4_-_Star_edge_{1,_4}: - 2 partition[4,4]
   - circuit_partition_node[1,4] + star_edge[1,4] >= -2
 Partition_4_-_Star_edge_{1,_4}: - 2 partition[1,4]
   - circuit_partition_node[4,4] + star_edge[1,4] >= -2
 Partition_4_-_Star_edge_{1,_3}: - 2 partition[3,4]
   - circuit_partition_node[1,4] + star_edge[1,3] >= -2
 Partition_4_-_Star_edge_{1,_3}: - 2 partition[1,4]
   - circuit_partition_node[3,4] + star_edge[1,3] >= -2
 Partition_4_-_Star_edge_{1,_2}: - 2 partition[2,4]
   - circuit_partition_node[1,4] + star_edge[1,2] >= -2
 Partition_4_-_Star_edge_{1,_2}: - 2 partition[1,4]
   - circuit_partition_node[2,4] + star_edge[1,2] >= -2
 Partition_4_-_Star_edge_{7,_8}: - 2 partition[8,4]
   - circuit_partition_node[7,4] + star_edge[7,8] >= -2
 Partition_4_-_Star_edge_{7,_8}: - 2 partition[7,4]
   - circuit_partition_node[8,4] + star_edge[7,8] >= -2
 Partition_4_-_Star_edge_{7,_9}: - 2 partition[9,4]
   - circuit_partition_node[7,4] + star_edge[7,9] >= -2
 Partition_4_-_Star_edge_{7,_9}: - 2 partition[7,4]
   - circuit_partition_node[9,4] + star_edge[7,9] >= -2
 Partition_4_-_Star_edge_{0,_9}: - 2 partition[9,4]
   - circuit_partition_node[0,4] + star_edge[0,9] >= -2
 Partition_4_-_Star_edge_{0,_9}: - 2 partition[0,4]
   - circuit_partition_node[9,4] + star_edge[0,9] >= -2
 Partition_4_-_Star_edge_{0,_8}: - 2 partition[8,4]
   - circuit_partition_node[0,4] + star_edge[0,8] >= -2
 Partition_4_-_Star_edge_{0,_8}: - 2 partition[0,4]
   - circuit_partition_node[8,4] + star_edge[0,8] >= -2
 Partition_4_-_Star_edge_{0,_7}: - 2 partition[7,4]
   - circuit_partition_node[0,4] + star_edge[0,7] >= -2
 Partition_4_-_Star_edge_{0,_7}: - 2 partition[0,4]
   - circuit_partition_node[7,4] + star_edge[0,7] >= -2
 Partition_4_-_Star_edge_{0,_6}: - 2 partition[6,4]
   - circuit_partition_node[0,4] + star_edge[0,6] >= -2
 Partition_4_-_Star_edge_{0,_6}: - 2 partition[0,4]
   - circuit_partition_node[6,4] + star_edge[0,6] >= -2
 Partition_4_-_Star_edge_{0,_5}: - 2 partition[5,4]
   - circuit_partition_node[0,4] + star_edge[0,5] >= -2
 Partition_4_-_Star_edge_{0,_5}: - 2 partition[0,4]
   - circuit_partition_node[5,4] + star_edge[0,5] >= -2
 Partition_4_-_Star_edge_{0,_4}: - 2 partition[4,4]
   - circuit_partition_node[0,4] + star_edge[0,4] >= -2
 Partition_4_-_Star_edge_{0,_4}: - 2 partition[0,4]
   - circuit_partition_node[4,4] + star_edge[0,4] >= -2
 Partition_4_-_Star_edge_{0,_3}: - 2 partition[3,4]
   - circuit_partition_node[0,4] + star_edge[0,3] >= -2
 Partition_4_-_Star_edge_{0,_3}: - 2 partition[0,4]
   - circuit_partition_node[3,4] + star_edge[0,3] >= -2
 Partition_4_-_Star_edge_{0,_2}: - 2 partition[2,4]
   - circuit_partition_node[0,4] + star_edge[0,2] >= -2
 Partition_4_-_Star_edge_{0,_2}: - 2 partition[0,4]
   - circuit_partition_node[2,4] + star_edge[0,2] >= -2
 Partition_4_-_Star_edge_{0,_1}: - 2 partition[1,4]
   - circuit_partition_node[0,4] + star_edge[0,1] >= -2
 Partition_4_-_Star_edge_{0,_1}: - 2 partition[0,4]
   - circuit_partition_node[1,4] + star_edge[0,1] >= -2
 Partition_4_-_Star_edge_{8,_9}: - 2 partition[9,4]
   - circuit_partition_node[8,4] + star_edge[8,9] >= -2
 Partition_4_-_Star_edge_{8,_9}: - 2 partition[8,4]
   - circuit_partition_node[9,4] + star_edge[8,9] >= -2
 Partition_5_-_Star_edge_{4,_5}: - 2 partition[5,5]
   - circuit_partition_node[4,5] + star_edge[4,5] >= -2
 Partition_5_-_Star_edge_{4,_5}: - 2 partition[4,5]
   - circuit_partition_node[5,5] + star_edge[4,5] >= -2
 Partition_5_-_Star_edge_{4,_6}: - 2 partition[6,5]
   - circuit_partition_node[4,5] + star_edge[4,6] >= -2
 Partition_5_-_Star_edge_{4,_6}: - 2 partition[4,5]
   - circuit_partition_node[6,5] + star_edge[4,6] >= -2
 Partition_5_-_Star_edge_{4,_7}: - 2 partition[7,5]
   - circuit_partition_node[4,5] + star_edge[4,7] >= -2
 Partition_5_-_Star_edge_{4,_7}: - 2 partition[4,5]
   - circuit_partition_node[7,5] + star_edge[4,7] >= -2
 Partition_5_-_Star_edge_{4,_8}: - 2 partition[8,5]
   - circuit_partition_node[4,5] + star_edge[4,8] >= -2
 Partition_5_-_Star_edge_{4,_8}: - 2 partition[4,5]
   - circuit_partition_node[8,5] + star_edge[4,8] >= -2
 Partition_5_-_Star_edge_{4,_9}: - 2 partition[9,5]
   - circuit_partition_node[4,5] + star_edge[4,9] >= -2
 Partition_5_-_Star_edge_{4,_9}: - 2 partition[4,5]
   - circuit_partition_node[9,5] + star_edge[4,9] >= -2
 Partition_5_-_Star_edge_{3,_9}: - 2 partition[9,5]
   - circuit_partition_node[3,5] + star_edge[3,9] >= -2
 Partition_5_-_Star_edge_{3,_9}: - 2 partition[3,5]
   - circuit_partition_node[9,5] + star_edge[3,9] >= -2
 Partition_5_-_Star_edge_{3,_8}: - 2 partition[8,5]
   - circuit_partition_node[3,5] + star_edge[3,8] >= -2
 Partition_5_-_Star_edge_{3,_8}: - 2 partition[3,5]
   - circuit_partition_node[8,5] + star_edge[3,8] >= -2
 Partition_5_-_Star_edge_{3,_7}: - 2 partition[7,5]
   - circuit_partition_node[3,5] + star_edge[3,7] >= -2
 Partition_5_-_Star_edge_{3,_7}: - 2 partition[3,5]
   - circuit_partition_node[7,5] + star_edge[3,7] >= -2
 Partition_5_-_Star_edge_{3,_6}: - 2 partition[6,5]
   - circuit_partition_node[3,5] + star_edge[3,6] >= -2
 Partition_5_-_Star_edge_{3,_6}: - 2 partition[3,5]
   - circuit_partition_node[6,5] + star_edge[3,6] >= -2
 Partition_5_-_Star_edge_{3,_5}: - 2 partition[5,5]
   - circuit_partition_node[3,5] + star_edge[3,5] >= -2
 Partition_5_-_Star_edge_{3,_5}: - 2 partition[3,5]
   - circuit_partition_node[5,5] + star_edge[3,5] >= -2
 Partition_5_-_Star_edge_{3,_4}: - 2 partition[4,5]
   - circuit_partition_node[3,5] + star_edge[3,4] >= -2
 Partition_5_-_Star_edge_{3,_4}: - 2 partition[3,5]
   - circuit_partition_node[4,5] + star_edge[3,4] >= -2
 Partition_5_-_Star_edge_{5,_6}: - 2 partition[6,5]
   - circuit_partition_node[5,5] + star_edge[5,6] >= -2
 Partition_5_-_Star_edge_{5,_6}: - 2 partition[5,5]
   - circuit_partition_node[6,5] + star_edge[5,6] >= -2
 Partition_5_-_Star_edge_{5,_7}: - 2 partition[7,5]
   - circuit_partition_node[5,5] + star_edge[5,7] >= -2
 Partition_5_-_Star_edge_{5,_7}: - 2 partition[5,5]
   - circuit_partition_node[7,5] + star_edge[5,7] >= -2
 Partition_5_-_Star_edge_{5,_8}: - 2 partition[8,5]
   - circuit_partition_node[5,5] + star_edge[5,8] >= -2
 Partition_5_-_Star_edge_{5,_8}: - 2 partition[5,5]
   - circuit_partition_node[8,5] + star_edge[5,8] >= -2
 Partition_5_-_Star_edge_{5,_9}: - 2 partition[9,5]
   - circuit_partition_node[5,5] + star_edge[5,9] >= -2
 Partition_5_-_Star_edge_{5,_9}: - 2 partition[5,5]
   - circuit_partition_node[9,5] + star_edge[5,9] >= -2
 Partition_5_-_Star_edge_{2,_9}: - 2 partition[9,5]
   - circuit_partition_node[2,5] + star_edge[2,9] >= -2
 Partition_5_-_Star_edge_{2,_9}: - 2 partition[2,5]
   - circuit_partition_node[9,5] + star_edge[2,9] >= -2
 Partition_5_-_Star_edge_{2,_8}: - 2 partition[8,5]
   - circuit_partition_node[2,5] + star_edge[2,8] >= -2
 Partition_5_-_Star_edge_{2,_8}: - 2 partition[2,5]
   - circuit_partition_node[8,5] + star_edge[2,8] >= -2
 Partition_5_-_Star_edge_{2,_7}: - 2 partition[7,5]
   - circuit_partition_node[2,5] + star_edge[2,7] >= -2
 Partition_5_-_Star_edge_{2,_7}: - 2 partition[2,5]
   - circuit_partition_node[7,5] + star_edge[2,7] >= -2
 Partition_5_-_Star_edge_{2,_6}: - 2 partition[6,5]
   - circuit_partition_node[2,5] + star_edge[2,6] >= -2
 Partition_5_-_Star_edge_{2,_6}: - 2 partition[2,5]
   - circuit_partition_node[6,5] + star_edge[2,6] >= -2
 Partition_5_-_Star_edge_{2,_5}: - 2 partition[5,5]
   - circuit_partition_node[2,5] + star_edge[2,5] >= -2
 Partition_5_-_Star_edge_{2,_5}: - 2 partition[2,5]
   - circuit_partition_node[5,5] + star_edge[2,5] >= -2
 Partition_5_-_Star_edge_{2,_4}: - 2 partition[4,5]
   - circuit_partition_node[2,5] + star_edge[2,4] >= -2
 Partition_5_-_Star_edge_{2,_4}: - 2 partition[2,5]
   - circuit_partition_node[4,5] + star_edge[2,4] >= -2
 Partition_5_-_Star_edge_{2,_3}: - 2 partition[3,5]
   - circuit_partition_node[2,5] + star_edge[2,3] >= -2
 Partition_5_-_Star_edge_{2,_3}: - 2 partition[2,5]
   - circuit_partition_node[3,5] + star_edge[2,3] >= -2
 Partition_5_-_Star_edge_{6,_7}: - 2 partition[7,5]
   - circuit_partition_node[6,5] + star_edge[6,7] >= -2
 Partition_5_-_Star_edge_{6,_7}: - 2 partition[6,5]
   - circuit_partition_node[7,5] + star_edge[6,7] >= -2
 Partition_5_-_Star_edge_{6,_8}: - 2 partition[8,5]
   - circuit_partition_node[6,5] + star_edge[6,8] >= -2
 Partition_5_-_Star_edge_{6,_8}: - 2 partition[6,5]
   - circuit_partition_node[8,5] + star_edge[6,8] >= -2
 Partition_5_-_Star_edge_{6,_9}: - 2 partition[9,5]
   - circuit_partition_node[6,5] + star_edge[6,9] >= -2
 Partition_5_-_Star_edge_{6,_9}: - 2 partition[6,5]
   - circuit_partition_node[9,5] + star_edge[6,9] >= -2
 Partition_5_-_Star_edge_{1,_9}: - 2 partition[9,5]
   - circuit_partition_node[1,5] + star_edge[1,9] >= -2
 Partition_5_-_Star_edge_{1,_9}: - 2 partition[1,5]
   - circuit_partition_node[9,5] + star_edge[1,9] >= -2
 Partition_5_-_Star_edge_{1,_8}: - 2 partition[8,5]
   - circuit_partition_node[1,5] + star_edge[1,8] >= -2
 Partition_5_-_Star_edge_{1,_8}: - 2 partition[1,5]
   - circuit_partition_node[8,5] + star_edge[1,8] >= -2
 Partition_5_-_Star_edge_{1,_7}: - 2 partition[7,5]
   - circuit_partition_node[1,5] + star_edge[1,7] >= -2
 Partition_5_-_Star_edge_{1,_7}: - 2 partition[1,5]
   - circuit_partition_node[7,5] + star_edge[1,7] >= -2
 Partition_5_-_Star_edge_{1,_6}: - 2 partition[6,5]
   - circuit_partition_node[1,5] + star_edge[1,6] >= -2
 Partition_5_-_Star_edge_{1,_6}: - 2 partition[1,5]
   - circuit_partition_node[6,5] + star_edge[1,6] >= -2
 Partition_5_-_Star_edge_{1,_5}: - 2 partition[5,5]
   - circuit_partition_node[1,5] + star_edge[1,5] >= -2
 Partition_5_-_Star_edge_{1,_5}: - 2 partition[1,5]
   - circuit_partition_node[5,5] + star_edge[1,5] >= -2
 Partition_5_-_Star_edge_{1,_4}: - 2 partition[4,5]
   - circuit_partition_node[1,5] + star_edge[1,4] >= -2
 Partition_5_-_Star_edge_{1,_4}: - 2 partition[1,5]
   - circuit_partition_node[4,5] + star_edge[1,4] >= -2
 Partition_5_-_Star_edge_{1,_3}: - 2 partition[3,5]
   - circuit_partition_node[1,5] + star_edge[1,3] >= -2
 Partition_5_-_Star_edge_{1,_3}: - 2 partition[1,5]
   - circuit_partition_node[3,5] + star_edge[1,3] >= -2
 Partition_5_-_Star_edge_{1,_2}: - 2 partition[2,5]
   - circuit_partition_node[1,5] + star_edge[1,2] >= -2
 Partition_5_-_Star_edge_{1,_2}: - 2 partition[1,5]
   - circuit_partition_node[2,5] + star_edge[1,2] >= -2
 Partition_5_-_Star_edge_{7,_8}: - 2 partition[8,5]
   - circuit_partition_node[7,5] + star_edge[7,8] >= -2
 Partition_5_-_Star_edge_{7,_8}: - 2 partition[7,5]
   - circuit_partition_node[8,5] + star_edge[7,8] >= -2
 Partition_5_-_Star_edge_{7,_9}: - 2 partition[9,5]
   - circuit_partition_node[7,5] + star_edge[7,9] >= -2
 Partition_5_-_Star_edge_{7,_9}: - 2 partition[7,5]
   - circuit_partition_node[9,5] + star_edge[7,9] >= -2
 Partition_5_-_Star_edge_{0,_9}: - 2 partition[9,5]
   - circuit_partition_node[0,5] + star_edge[0,9] >= -2
 Partition_5_-_Star_edge_{0,_9}: - 2 partition[0,5]
   - circuit_partition_node[9,5] + star_edge[0,9] >= -2
 Partition_5_-_Star_edge_{0,_8}: - 2 partition[8,5]
   - circuit_partition_node[0,5] + star_edge[0,8] >= -2
 Partition_5_-_Star_edge_{0,_8}: - 2 partition[0,5]
   - circuit_partition_node[8,5] + star_edge[0,8] >= -2
 Partition_5_-_Star_edge_{0,_7}: - 2 partition[7,5]
   - circuit_partition_node[0,5] + star_edge[0,7] >= -2
 Partition_5_-_Star_edge_{0,_7}: - 2 partition[0,5]
   - circuit_partition_node[7,5] + star_edge[0,7] >= -2
 Partition_5_-_Star_edge_{0,_6}: - 2 partition[6,5]
   - circuit_partition_node[0,5] + star_edge[0,6] >= -2
 Partition_5_-_Star_edge_{0,_6}: - 2 partition[0,5]
   - circuit_partition_node[6,5] + star_edge[0,6] >= -2
 Partition_5_-_Star_edge_{0,_5}: - 2 partition[5,5]
   - circuit_partition_node[0,5] + star_edge[0,5] >= -2
 Partition_5_-_Star_edge_{0,_5}: - 2 partition[0,5]
   - circuit_partition_node[5,5] + star_edge[0,5] >= -2
 Partition_5_-_Star_edge_{0,_4}: - 2 partition[4,5]
   - circuit_partition_node[0,5] + star_edge[0,4] >= -2
 Partition_5_-_Star_edge_{0,_4}: - 2 partition[0,5]
   - circuit_partition_node[4,5] + star_edge[0,4] >= -2
 Partition_5_-_Star_edge_{0,_3}: - 2 partition[3,5]
   - circuit_partition_node[0,5] + star_edge[0,3] >= -2
 Partition_5_-_Star_edge_{0,_3}: - 2 partition[0,5]
   - circuit_partition_node[3,5] + star_edge[0,3] >= -2
 Partition_5_-_Star_edge_{0,_2}: - 2 partition[2,5]
   - circuit_partition_node[0,5] + star_edge[0,2] >= -2
 Partition_5_-_Star_edge_{0,_2}: - 2 partition[0,5]
   - circuit_partition_node[2,5] + star_edge[0,2] >= -2
 Partition_5_-_Star_edge_{0,_1}: - 2 partition[1,5]
   - circuit_partition_node[0,5] + star_edge[0,1] >= -2
 Partition_5_-_Star_edge_{0,_1}: - 2 partition[0,5]
   - circuit_partition_node[1,5] + star_edge[0,1] >= -2
 Partition_5_-_Star_edge_{8,_9}: - 2 partition[9,5]
   - circuit_partition_node[8,5] + star_edge[8,9] >= -2
 Partition_5_-_Star_edge_{8,_9}: - 2 partition[8,5]
   - circuit_partition_node[9,5] + star_edge[8,9] >= -2
 Partition_6_-_Star_edge_{4,_5}: - 2 partition[5,6]
   - circuit_partition_node[4,6] + star_edge[4,5] >= -2
 Partition_6_-_Star_edge_{4,_5}: - 2 partition[4,6]
   - circuit_partition_node[5,6] + star_edge[4,5] >= -2
 Partition_6_-_Star_edge_{4,_6}: - 2 partition[6,6]
   - circuit_partition_node[4,6] + star_edge[4,6] >= -2
 Partition_6_-_Star_edge_{4,_6}: - 2 partition[4,6]
   - circuit_partition_node[6,6] + star_edge[4,6] >= -2
 Partition_6_-_Star_edge_{4,_7}: - 2 partition[7,6]
   - circuit_partition_node[4,6] + star_edge[4,7] >= -2
 Partition_6_-_Star_edge_{4,_7}: - 2 partition[4,6]
   - circuit_partition_node[7,6] + star_edge[4,7] >= -2
 Partition_6_-_Star_edge_{4,_8}: - 2 partition[8,6]
   - circuit_partition_node[4,6] + star_edge[4,8] >= -2
 Partition_6_-_Star_edge_{4,_8}: - 2 partition[4,6]
   - circuit_partition_node[8,6] + star_edge[4,8] >= -2
 Partition_6_-_Star_edge_{4,_9}: - 2 partition[9,6]
   - circuit_partition_node[4,6] + star_edge[4,9] >= -2
 Partition_6_-_Star_edge_{4,_9}: - 2 partition[4,6]
   - circuit_partition_node[9,6] + star_edge[4,9] >= -2
 Partition_6_-_Star_edge_{3,_9}: - 2 partition[9,6]
   - circuit_partition_node[3,6] + star_edge[3,9] >= -2
 Partition_6_-_Star_edge_{3,_9}: - 2 partition[3,6]
   - circuit_partition_node[9,6] + star_edge[3,9] >= -2
 Partition_6_-_Star_edge_{3,_8}: - 2 partition[8,6]
   - circuit_partition_node[3,6] + star_edge[3,8] >= -2
 Partition_6_-_Star_edge_{3,_8}: - 2 partition[3,6]
   - circuit_partition_node[8,6] + star_edge[3,8] >= -2
 Partition_6_-_Star_edge_{3,_7}: - 2 partition[7,6]
   - circuit_partition_node[3,6] + star_edge[3,7] >= -2
 Partition_6_-_Star_edge_{3,_7}: - 2 partition[3,6]
   - circuit_partition_node[7,6] + star_edge[3,7] >= -2
 Partition_6_-_Star_edge_{3,_6}: - 2 partition[6,6]
   - circuit_partition_node[3,6] + star_edge[3,6] >= -2
 Partition_6_-_Star_edge_{3,_6}: - 2 partition[3,6]
   - circuit_partition_node[6,6] + star_edge[3,6] >= -2
 Partition_6_-_Star_edge_{3,_5}: - 2 partition[5,6]
   - circuit_partition_node[3,6] + star_edge[3,5] >= -2
 Partition_6_-_Star_edge_{3,_5}: - 2 partition[3,6]
   - circuit_partition_node[5,6] + star_edge[3,5] >= -2
 Partition_6_-_Star_edge_{3,_4}: - 2 partition[4,6]
   - circuit_partition_node[3,6] + star_edge[3,4] >= -2
 Partition_6_-_Star_edge_{3,_4}: - 2 partition[3,6]
   - circuit_partition_node[4,6] + star_edge[3,4] >= -2
 Partition_6_-_Star_edge_{5,_6}: - 2 partition[6,6]
   - circuit_partition_node[5,6] + star_edge[5,6] >= -2
 Partition_6_-_Star_edge_{5,_6}: - 2 partition[5,6]
   - circuit_partition_node[6,6] + star_edge[5,6] >= -2
 Partition_6_-_Star_edge_{5,_7}: - 2 partition[7,6]
   - circuit_partition_node[5,6] + star_edge[5,7] >= -2
 Partition_6_-_Star_edge_{5,_7}: - 2 partition[5,6]
   - circuit_partition_node[7,6] + star_edge[5,7] >= -2
 Partition_6_-_Star_edge_{5,_8}: - 2 partition[8,6]
   - circuit_partition_node[5,6] + star_edge[5,8] >= -2
 Partition_6_-_Star_edge_{5,_8}: - 2 partition[5,6]
   - circuit_partition_node[8,6] + star_edge[5,8] >= -2
 Partition_6_-_Star_edge_{5,_9}: - 2 partition[9,6]
   - circuit_partition_node[5,6] + star_edge[5,9] >= -2
 Partition_6_-_Star_edge_{5,_9}: - 2 partition[5,6]
   - circuit_partition_node[9,6] + star_edge[5,9] >= -2
 Partition_6_-_Star_edge_{2,_9}: - 2 partition[9,6]
   - circuit_partition_node[2,6] + star_edge[2,9] >= -2
 Partition_6_-_Star_edge_{2,_9}: - 2 partition[2,6]
   - circuit_partition_node[9,6] + star_edge[2,9] >= -2
 Partition_6_-_Star_edge_{2,_8}: - 2 partition[8,6]
   - circuit_partition_node[2,6] + star_edge[2,8] >= -2
 Partition_6_-_Star_edge_{2,_8}: - 2 partition[2,6]
   - circuit_partition_node[8,6] + star_edge[2,8] >= -2
 Partition_6_-_Star_edge_{2,_7}: - 2 partition[7,6]
   - circuit_partition_node[2,6] + star_edge[2,7] >= -2
 Partition_6_-_Star_edge_{2,_7}: - 2 partition[2,6]
   - circuit_partition_node[7,6] + star_edge[2,7] >= -2
 Partition_6_-_Star_edge_{2,_6}: - 2 partition[6,6]
   - circuit_partition_node[2,6] + star_edge[2,6] >= -2
 Partition_6_-_Star_edge_{2,_6}: - 2 partition[2,6]
   - circuit_partition_node[6,6] + star_edge[2,6] >= -2
 Partition_6_-_Star_edge_{2,_5}: - 2 partition[5,6]
   - circuit_partition_node[2,6] + star_edge[2,5] >= -2
 Partition_6_-_Star_edge_{2,_5}: - 2 partition[2,6]
   - circuit_partition_node[5,6] + star_edge[2,5] >= -2
 Partition_6_-_Star_edge_{2,_4}: - 2 partition[4,6]
   - circuit_partition_node[2,6] + star_edge[2,4] >= -2
 Partition_6_-_Star_edge_{2,_4}: - 2 partition[2,6]
   - circuit_partition_node[4,6] + star_edge[2,4] >= -2
 Partition_6_-_Star_edge_{2,_3}: - 2 partition[3,6]
   - circuit_partition_node[2,6] + star_edge[2,3] >= -2
 Partition_6_-_Star_edge_{2,_3}: - 2 partition[2,6]
   - circuit_partition_node[3,6] + star_edge[2,3] >= -2
 Partition_6_-_Star_edge_{6,_7}: - 2 partition[7,6]
   - circuit_partition_node[6,6] + star_edge[6,7] >= -2
 Partition_6_-_Star_edge_{6,_7}: - 2 partition[6,6]
   - circuit_partition_node[7,6] + star_edge[6,7] >= -2
 Partition_6_-_Star_edge_{6,_8}: - 2 partition[8,6]
   - circuit_partition_node[6,6] + star_edge[6,8] >= -2
 Partition_6_-_Star_edge_{6,_8}: - 2 partition[6,6]
   - circuit_partition_node[8,6] + star_edge[6,8] >= -2
 Partition_6_-_Star_edge_{6,_9}: - 2 partition[9,6]
   - circuit_partition_node[6,6] + star_edge[6,9] >= -2
 Partition_6_-_Star_edge_{6,_9}: - 2 partition[6,6]
   - circuit_partition_node[9,6] + star_edge[6,9] >= -2
 Partition_6_-_Star_edge_{1,_9}: - 2 partition[9,6]
   - circuit_partition_node[1,6] + star_edge[1,9] >= -2
 Partition_6_-_Star_edge_{1,_9}: - 2 partition[1,6]
   - circuit_partition_node[9,6] + star_edge[1,9] >= -2
 Partition_6_-_Star_edge_{1,_8}: - 2 partition[8,6]
   - circuit_partition_node[1,6] + star_edge[1,8] >= -2
 Partition_6_-_Star_edge_{1,_8}: - 2 partition[1,6]
   - circuit_partition_node[8,6] + star_edge[1,8] >= -2
 Partition_6_-_Star_edge_{1,_7}: - 2 partition[7,6]
   - circuit_partition_node[1,6] + star_edge[1,7] >= -2
 Partition_6_-_Star_edge_{1,_7}: - 2 partition[1,6]
   - circuit_partition_node[7,6] + star_edge[1,7] >= -2
 Partition_6_-_Star_edge_{1,_6}: - 2 partition[6,6]
   - circuit_partition_node[1,6] + star_edge[1,6] >= -2
 Partition_6_-_Star_edge_{1,_6}: - 2 partition[1,6]
   - circuit_partition_node[6,6] + star_edge[1,6] >= -2
 Partition_6_-_Star_edge_{1,_5}: - 2 partition[5,6]
   - circuit_partition_node[1,6] + star_edge[1,5] >= -2
 Partition_6_-_Star_edge_{1,_5}: - 2 partition[1,6]
   - circuit_partition_node[5,6] + star_edge[1,5] >= -2
 Partition_6_-_Star_edge_{1,_4}: - 2 partition[4,6]
   - circuit_partition_node[1,6] + star_edge[1,4] >= -2
 Partition_6_-_Star_edge_{1,_4}: - 2 partition[1,6]
   - circuit_partition_node[4,6] + star_edge[1,4] >= -2
 Partition_6_-_Star_edge_{1,_3}: - 2 partition[3,6]
   - circuit_partition_node[1,6] + star_edge[1,3] >= -2
 Partition_6_-_Star_edge_{1,_3}: - 2 partition[1,6]
   - circuit_partition_node[3,6] + star_edge[1,3] >= -2
 Partition_6_-_Star_edge_{1,_2}: - 2 partition[2,6]
   - circuit_partition_node[1,6] + star_edge[1,2] >= -2
 Partition_6_-_Star_edge_{1,_2}: - 2 partition[1,6]
   - circuit_partition_node[2,6] + star_edge[1,2] >= -2
 Partition_6_-_Star_edge_{7,_8}: - 2 partition[8,6]
   - circuit_partition_node[7,6] + star_edge[7,8] >= -2
 Partition_6_-_Star_edge_{7,_8}: - 2 partition[7,6]
   - circuit_partition_node[8,6] + star_edge[7,8] >= -2
 Partition_6_-_Star_edge_{7,_9}: - 2 partition[9,6]
   - circuit_partition_node[7,6] + star_edge[7,9] >= -2
 Partition_6_-_Star_edge_{7,_9}: - 2 partition[7,6]
   - circuit_partition_node[9,6] + star_edge[7,9] >= -2
 Partition_6_-_Star_edge_{0,_9}: - 2 partition[9,6]
   - circuit_partition_node[0,6] + star_edge[0,9] >= -2
 Partition_6_-_Star_edge_{0,_9}: - 2 partition[0,6]
   - circuit_partition_node[9,6] + star_edge[0,9] >= -2
 Partition_6_-_Star_edge_{0,_8}: - 2 partition[8,6]
   - circuit_partition_node[0,6] + star_edge[0,8] >= -2
 Partition_6_-_Star_edge_{0,_8}: - 2 partition[0,6]
   - circuit_partition_node[8,6] + star_edge[0,8] >= -2
 Partition_6_-_Star_edge_{0,_7}: - 2 partition[7,6]
   - circuit_partition_node[0,6] + star_edge[0,7] >= -2
 Partition_6_-_Star_edge_{0,_7}: - 2 partition[0,6]
   - circuit_partition_node[7,6] + star_edge[0,7] >= -2
 Partition_6_-_Star_edge_{0,_6}: - 2 partition[6,6]
   - circuit_partition_node[0,6] + star_edge[0,6] >= -2
 Partition_6_-_Star_edge_{0,_6}: - 2 partition[0,6]
   - circuit_partition_node[6,6] + star_edge[0,6] >= -2
 Partition_6_-_Star_edge_{0,_5}: - 2 partition[5,6]
   - circuit_partition_node[0,6] + star_edge[0,5] >= -2
 Partition_6_-_Star_edge_{0,_5}: - 2 partition[0,6]
   - circuit_partition_node[5,6] + star_edge[0,5] >= -2
 Partition_6_-_Star_edge_{0,_4}: - 2 partition[4,6]
   - circuit_partition_node[0,6] + star_edge[0,4] >= -2
 Partition_6_-_Star_edge_{0,_4}: - 2 partition[0,6]
   - circuit_partition_node[4,6] + star_edge[0,4] >= -2
 Partition_6_-_Star_edge_{0,_3}: - 2 partition[3,6]
   - circuit_partition_node[0,6] + star_edge[0,3] >= -2
 Partition_6_-_Star_edge_{0,_3}: - 2 partition[0,6]
   - circuit_partition_node[3,6] + star_edge[0,3] >= -2
 Partition_6_-_Star_edge_{0,_2}: - 2 partition[2,6]
   - circuit_partition_node[0,6] + star_edge[0,2] >= -2
 Partition_6_-_Star_edge_{0,_2}: - 2 partition[0,6]
   - circuit_partition_node[2,6] + star_edge[0,2] >= -2
 Partition_6_-_Star_edge_{0,_1}: - 2 partition[1,6]
   - circuit_partition_node[0,6] + star_edge[0,1] >= -2
 Partition_6_-_Star_edge_{0,_1}: - 2 partition[0,6]
   - circuit_partition_node[1,6] + star_edge[0,1] >= -2
 Partition_6_-_Star_edge_{8,_9}: - 2 partition[9,6]
   - circuit_partition_node[8,6] + star_edge[8,9] >= -2
 Partition_6_-_Star_edge_{8,_9}: - 2 partition[8,6]
   - circuit_partition_node[9,6] + star_edge[8,9] >= -2
 Partition_7_-_Star_edge_{4,_5}: - 2 partition[5,7]
   - circuit_partition_node[4,7] + star_edge[4,5] >= -2
 Partition_7_-_Star_edge_{4,_5}: - 2 partition[4,7]
   - circuit_partition_node[5,7] + star_edge[4,5] >= -2
 Partition_7_-_Star_edge_{4,_6}: - 2 partition[6,7]
   - circuit_partition_node[4,7] + star_edge[4,6] >= -2
 Partition_7_-_Star_edge_{4,_6}: - 2 partition[4,7]
   - circuit_partition_node[6,7] + star_edge[4,6] >= -2
 Partition_7_-_Star_edge_{4,_7}: - 2 partition[7,7]
   - circuit_partition_node[4,7] + star_edge[4,7] >= -2
 Partition_7_-_Star_edge_{4,_7}: - 2 partition[4,7]
   - circuit_partition_node[7,7] + star_edge[4,7] >= -2
 Partition_7_-_Star_edge_{4,_8}: - 2 partition[8,7]
   - circuit_partition_node[4,7] + star_edge[4,8] >= -2
 Partition_7_-_Star_edge_{4,_8}: - 2 partition[4,7]
   - circuit_partition_node[8,7] + star_edge[4,8] >= -2
 Partition_7_-_Star_edge_{4,_9}: - 2 partition[9,7]
   - circuit_partition_node[4,7] + star_edge[4,9] >= -2
 Partition_7_-_Star_edge_{4,_9}: - 2 partition[4,7]
   - circuit_partition_node[9,7] + star_edge[4,9] >= -2
 Partition_7_-_Star_edge_{3,_9}: - 2 partition[9,7]
   - circuit_partition_node[3,7] + star_edge[3,9] >= -2
 Partition_7_-_Star_edge_{3,_9}: - 2 partition[3,7]
   - circuit_partition_node[9,7] + star_edge[3,9] >= -2
 Partition_7_-_Star_edge_{3,_8}: - 2 partition[8,7]
   - circuit_partition_node[3,7] + star_edge[3,8] >= -2
 Partition_7_-_Star_edge_{3,_8}: - 2 partition[3,7]
   - circuit_partition_node[8,7] + star_edge[3,8] >= -2
 Partition_7_-_Star_edge_{3,_7}: - 2 partition[7,7]
   - circuit_partition_node[3,7] + star_edge[3,7] >= -2
 Partition_7_-_Star_edge_{3,_7}: - 2 partition[3,7]
   - circuit_partition_node[7,7] + star_edge[3,7] >= -2
 Partition_7_-_Star_edge_{3,_6}: - 2 partition[6,7]
   - circuit_partition_node[3,7] + star_edge[3,6] >= -2
 Partition_7_-_Star_edge_{3,_6}: - 2 partition[3,7]
   - circuit_partition_node[6,7] + star_edge[3,6] >= -2
 Partition_7_-_Star_edge_{3,_5}: - 2 partition[5,7]
   - circuit_partition_node[3,7] + star_edge[3,5] >= -2
 Partition_7_-_Star_edge_{3,_5}: - 2 partition[3,7]
   - circuit_partition_node[5,7] + star_edge[3,5] >= -2
 Partition_7_-_Star_edge_{3,_4}: - 2 partition[4,7]
   - circuit_partition_node[3,7] + star_edge[3,4] >= -2
 Partition_7_-_Star_edge_{3,_4}: - 2 partition[3,7]
   - circuit_partition_node[4,7] + star_edge[3,4] >= -2
 Partition_7_-_Star_edge_{5,_6}: - 2 partition[6,7]
   - circuit_partition_node[5,7] + star_edge[5,6] >= -2
 Partition_7_-_Star_edge_{5,_6}: - 2 partition[5,7]
   - circuit_partition_node[6,7] + star_edge[5,6] >= -2
 Partition_7_-_Star_edge_{5,_7}: - 2 partition[7,7]
   - circuit_partition_node[5,7] + star_edge[5,7] >= -2
 Partition_7_-_Star_edge_{5,_7}: - 2 partition[5,7]
   - circuit_partition_node[7,7] + star_edge[5,7] >= -2
 Partition_7_-_Star_edge_{5,_8}: - 2 partition[8,7]
   - circuit_partition_node[5,7] + star_edge[5,8] >= -2
 Partition_7_-_Star_edge_{5,_8}: - 2 partition[5,7]
   - circuit_partition_node[8,7] + star_edge[5,8] >= -2
 Partition_7_-_Star_edge_{5,_9}: - 2 partition[9,7]
   - circuit_partition_node[5,7] + star_edge[5,9] >= -2
 Partition_7_-_Star_edge_{5,_9}: - 2 partition[5,7]
   - circuit_partition_node[9,7] + star_edge[5,9] >= -2
 Partition_7_-_Star_edge_{2,_9}: - 2 partition[9,7]
   - circuit_partition_node[2,7] + star_edge[2,9] >= -2
 Partition_7_-_Star_edge_{2,_9}: - 2 partition[2,7]
   - circuit_partition_node[9,7] + star_edge[2,9] >= -2
 Partition_7_-_Star_edge_{2,_8}: - 2 partition[8,7]
   - circuit_partition_node[2,7] + star_edge[2,8] >= -2
 Partition_7_-_Star_edge_{2,_8}: - 2 partition[2,7]
   - circuit_partition_node[8,7] + star_edge[2,8] >= -2
 Partition_7_-_Star_edge_{2,_7}: - 2 partition[7,7]
   - circuit_partition_node[2,7] + star_edge[2,7] >= -2
 Partition_7_-_Star_edge_{2,_7}: - 2 partition[2,7]
   - circuit_partition_node[7,7] + star_edge[2,7] >= -2
 Partition_7_-_Star_edge_{2,_6}: - 2 partition[6,7]
   - circuit_partition_node[2,7] + star_edge[2,6] >= -2
 Partition_7_-_Star_edge_{2,_6}: - 2 partition[2,7]
   - circuit_partition_node[6,7] + star_edge[2,6] >= -2
 Partition_7_-_Star_edge_{2,_5}: - 2 partition[5,7]
   - circuit_partition_node[2,7] + star_edge[2,5] >= -2
 Partition_7_-_Star_edge_{2,_5}: - 2 partition[2,7]
   - circuit_partition_node[5,7] + star_edge[2,5] >= -2
 Partition_7_-_Star_edge_{2,_4}: - 2 partition[4,7]
   - circuit_partition_node[2,7] + star_edge[2,4] >= -2
 Partition_7_-_Star_edge_{2,_4}: - 2 partition[2,7]
   - circuit_partition_node[4,7] + star_edge[2,4] >= -2
 Partition_7_-_Star_edge_{2,_3}: - 2 partition[3,7]
   - circuit_partition_node[2,7] + star_edge[2,3] >= -2
 Partition_7_-_Star_edge_{2,_3}: - 2 partition[2,7]
   - circuit_partition_node[3,7] + star_edge[2,3] >= -2
 Partition_7_-_Star_edge_{6,_7}: - 2 partition[7,7]
   - circuit_partition_node[6,7] + star_edge[6,7] >= -2
 Partition_7_-_Star_edge_{6,_7}: - 2 partition[6,7]
   - circuit_partition_node[7,7] + star_edge[6,7] >= -2
 Partition_7_-_Star_edge_{6,_8}: - 2 partition[8,7]
   - circuit_partition_node[6,7] + star_edge[6,8] >= -2
 Partition_7_-_Star_edge_{6,_8}: - 2 partition[6,7]
   - circuit_partition_node[8,7] + star_edge[6,8] >= -2
 Partition_7_-_Star_edge_{6,_9}: - 2 partition[9,7]
   - circuit_partition_node[6,7] + star_edge[6,9] >= -2
 Partition_7_-_Star_edge_{6,_9}: - 2 partition[6,7]
   - circuit_partition_node[9,7] + star_edge[6,9] >= -2
 Partition_7_-_Star_edge_{1,_9}: - 2 partition[9,7]
   - circuit_partition_node[1,7] + star_edge[1,9] >= -2
 Partition_7_-_Star_edge_{1,_9}: - 2 partition[1,7]
   - circuit_partition_node[9,7] + star_edge[1,9] >= -2
 Partition_7_-_Star_edge_{1,_8}: - 2 partition[8,7]
   - circuit_partition_node[1,7] + star_edge[1,8] >= -2
 Partition_7_-_Star_edge_{1,_8}: - 2 partition[1,7]
   - circuit_partition_node[8,7] + star_edge[1,8] >= -2
 Partition_7_-_Star_edge_{1,_7}: - 2 partition[7,7]
   - circuit_partition_node[1,7] + star_edge[1,7] >= -2
 Partition_7_-_Star_edge_{1,_7}: - 2 partition[1,7]
   - circuit_partition_node[7,7] + star_edge[1,7] >= -2
 Partition_7_-_Star_edge_{1,_6}: - 2 partition[6,7]
   - circuit_partition_node[1,7] + star_edge[1,6] >= -2
 Partition_7_-_Star_edge_{1,_6}: - 2 partition[1,7]
   - circuit_partition_node[6,7] + star_edge[1,6] >= -2
 Partition_7_-_Star_edge_{1,_5}: - 2 partition[5,7]
   - circuit_partition_node[1,7] + star_edge[1,5] >= -2
 Partition_7_-_Star_edge_{1,_5}: - 2 partition[1,7]
   - circuit_partition_node[5,7] + star_edge[1,5] >= -2
 Partition_7_-_Star_edge_{1,_4}: - 2 partition[4,7]
   - circuit_partition_node[1,7] + star_edge[1,4] >= -2
 Partition_7_-_Star_edge_{1,_4}: - 2 partition[1,7]
   - circuit_partition_node[4,7] + star_edge[1,4] >= -2
 Partition_7_-_Star_edge_{1,_3}: - 2 partition[3,7]
   - circuit_partition_node[1,7] + star_edge[1,3] >= -2
 Partition_7_-_Star_edge_{1,_3}: - 2 partition[1,7]
   - circuit_partition_node[3,7] + star_edge[1,3] >= -2
 Partition_7_-_Star_edge_{1,_2}: - 2 partition[2,7]
   - circuit_partition_node[1,7] + star_edge[1,2] >= -2
 Partition_7_-_Star_edge_{1,_2}: - 2 partition[1,7]
   - circuit_partition_node[2,7] + star_edge[1,2] >= -2
 Partition_7_-_Star_edge_{7,_8}: - 2 partition[8,7]
   - circuit_partition_node[7,7] + star_edge[7,8] >= -2
 Partition_7_-_Star_edge_{7,_8}: - 2 partition[7,7]
   - circuit_partition_node[8,7] + star_edge[7,8] >= -2
 Partition_7_-_Star_edge_{7,_9}: - 2 partition[9,7]
   - circuit_partition_node[7,7] + star_edge[7,9] >= -2
 Partition_7_-_Star_edge_{7,_9}: - 2 partition[7,7]
   - circuit_partition_node[9,7] + star_edge[7,9] >= -2
 Partition_7_-_Star_edge_{0,_9}: - 2 partition[9,7]
   - circuit_partition_node[0,7] + star_edge[0,9] >= -2
 Partition_7_-_Star_edge_{0,_9}: - 2 partition[0,7]
   - circuit_partition_node[9,7] + star_edge[0,9] >= -2
 Partition_7_-_Star_edge_{0,_8}: - 2 partition[8,7]
   - circuit_partition_node[0,7] + star_edge[0,8] >= -2
 Partition_7_-_Star_edge_{0,_8}: - 2 partition[0,7]
   - circuit_partition_node[8,7] + star_edge[0,8] >= -2
 Partition_7_-_Star_edge_{0,_7}: - 2 partition[7,7]
   - circuit_partition_node[0,7] + star_edge[0,7] >= -2
 Partition_7_-_Star_edge_{0,_7}: - 2 partition[0,7]
   - circuit_partition_node[7,7] + star_edge[0,7] >= -2
 Partition_7_-_Star_edge_{0,_6}: - 2 partition[6,7]
   - circuit_partition_node[0,7] + star_edge[0,6] >= -2
 Partition_7_-_Star_edge_{0,_6}: - 2 partition[0,7]
   - circuit_partition_node[6,7] + star_edge[0,6] >= -2
 Partition_7_-_Star_edge_{0,_5}: - 2 partition[5,7]
   - circuit_partition_node[0,7] + star_edge[0,5] >= -2
 Partition_7_-_Star_edge_{0,_5}: - 2 partition[0,7]
   - circuit_partition_node[5,7] + star_edge[0,5] >= -2
 Partition_7_-_Star_edge_{0,_4}: - 2 partition[4,7]
   - circuit_partition_node[0,7] + star_edge[0,4] >= -2
 Partition_7_-_Star_edge_{0,_4}: - 2 partition[0,7]
   - circuit_partition_node[4,7] + star_edge[0,4] >= -2
 Partition_7_-_Star_edge_{0,_3}: - 2 partition[3,7]
   - circuit_partition_node[0,7] + star_edge[0,3] >= -2
 Partition_7_-_Star_edge_{0,_3}: - 2 partition[0,7]
   - circuit_partition_node[3,7] + star_edge[0,3] >= -2
 Partition_7_-_Star_edge_{0,_2}: - 2 partition[2,7]
   - circuit_partition_node[0,7] + star_edge[0,2] >= -2
 Partition_7_-_Star_edge_{0,_2}: - 2 partition[0,7]
   - circuit_partition_node[2,7] + star_edge[0,2] >= -2
 Partition_7_-_Star_edge_{0,_1}: - 2 partition[1,7]
   - circuit_partition_node[0,7] + star_edge[0,1] >= -2
 Partition_7_-_Star_edge_{0,_1}: - 2 partition[0,7]
   - circuit_partition_node[1,7] + star_edge[0,1] >= -2
 Partition_7_-_Star_edge_{8,_9}: - 2 partition[9,7]
   - circuit_partition_node[8,7] + star_edge[8,9] >= -2
 Partition_7_-_Star_edge_{8,_9}: - 2 partition[8,7]
   - circuit_partition_node[9,7] + star_edge[8,9] >= -2
 Partition_8_-_Star_edge_{4,_5}: - 2 partition[5,8]
   - circuit_partition_node[4,8] + star_edge[4,5] >= -2
 Partition_8_-_Star_edge_{4,_5}: - 2 partition[4,8]
   - circuit_partition_node[5,8] + star_edge[4,5] >= -2
 Partition_8_-_Star_edge_{4,_6}: - 2 partition[6,8]
   - circuit_partition_node[4,8] + star_edge[4,6] >= -2
 Partition_8_-_Star_edge_{4,_6}: - 2 partition[4,8]
   - circuit_partition_node[6,8] + star_edge[4,6] >= -2
 Partition_8_-_Star_edge_{4,_7}: - 2 partition[7,8]
   - circuit_partition_node[4,8] + star_edge[4,7] >= -2
 Partition_8_-_Star_edge_{4,_7}: - 2 partition[4,8]
   - circuit_partition_node[7,8] + star_edge[4,7] >= -2
 Partition_8_-_Star_edge_{4,_8}: - 2 partition[8,8]
   - circuit_partition_node[4,8] + star_edge[4,8] >= -2
 Partition_8_-_Star_edge_{4,_8}: - 2 partition[4,8]
   - circuit_partition_node[8,8] + star_edge[4,8] >= -2
 Partition_8_-_Star_edge_{4,_9}: - 2 partition[9,8]
   - circuit_partition_node[4,8] + star_edge[4,9] >= -2
 Partition_8_-_Star_edge_{4,_9}: - 2 partition[4,8]
   - circuit_partition_node[9,8] + star_edge[4,9] >= -2
 Partition_8_-_Star_edge_{3,_9}: - 2 partition[9,8]
   - circuit_partition_node[3,8] + star_edge[3,9] >= -2
 Partition_8_-_Star_edge_{3,_9}: - 2 partition[3,8]
   - circuit_partition_node[9,8] + star_edge[3,9] >= -2
 Partition_8_-_Star_edge_{3,_8}: - 2 partition[8,8]
   - circuit_partition_node[3,8] + star_edge[3,8] >= -2
 Partition_8_-_Star_edge_{3,_8}: - 2 partition[3,8]
   - circuit_partition_node[8,8] + star_edge[3,8] >= -2
 Partition_8_-_Star_edge_{3,_7}: - 2 partition[7,8]
   - circuit_partition_node[3,8] + star_edge[3,7] >= -2
 Partition_8_-_Star_edge_{3,_7}: - 2 partition[3,8]
   - circuit_partition_node[7,8] + star_edge[3,7] >= -2
 Partition_8_-_Star_edge_{3,_6}: - 2 partition[6,8]
   - circuit_partition_node[3,8] + star_edge[3,6] >= -2
 Partition_8_-_Star_edge_{3,_6}: - 2 partition[3,8]
   - circuit_partition_node[6,8] + star_edge[3,6] >= -2
 Partition_8_-_Star_edge_{3,_5}: - 2 partition[5,8]
   - circuit_partition_node[3,8] + star_edge[3,5] >= -2
 Partition_8_-_Star_edge_{3,_5}: - 2 partition[3,8]
   - circuit_partition_node[5,8] + star_edge[3,5] >= -2
 Partition_8_-_Star_edge_{3,_4}: - 2 partition[4,8]
   - circuit_partition_node[3,8] + star_edge[3,4] >= -2
 Partition_8_-_Star_edge_{3,_4}: - 2 partition[3,8]
   - circuit_partition_node[4,8] + star_edge[3,4] >= -2
 Partition_8_-_Star_edge_{5,_6}: - 2 partition[6,8]
   - circuit_partition_node[5,8] + star_edge[5,6] >= -2
 Partition_8_-_Star_edge_{5,_6}: - 2 partition[5,8]
   - circuit_partition_node[6,8] + star_edge[5,6] >= -2
 Partition_8_-_Star_edge_{5,_7}: - 2 partition[7,8]
   - circuit_partition_node[5,8] + star_edge[5,7] >= -2
 Partition_8_-_Star_edge_{5,_7}: - 2 partition[5,8]
   - circuit_partition_node[7,8] + star_edge[5,7] >= -2
 Partition_8_-_Star_edge_{5,_8}: - 2 partition[8,8]
   - circuit_partition_node[5,8] + star_edge[5,8] >= -2
 Partition_8_-_Star_edge_{5,_8}: - 2 partition[5,8]
   - circuit_partition_node[8,8] + star_edge[5,8] >= -2
 Partition_8_-_Star_edge_{5,_9}: - 2 partition[9,8]
   - circuit_partition_node[5,8] + star_edge[5,9] >= -2
 Partition_8_-_Star_edge_{5,_9}: - 2 partition[5,8]
   - circuit_partition_node[9,8] + star_edge[5,9] >= -2
 Partition_8_-_Star_edge_{2,_9}: - 2 partition[9,8]
   - circuit_partition_node[2,8] + star_edge[2,9] >= -2
 Partition_8_-_Star_edge_{2,_9}: - 2 partition[2,8]
   - circuit_partition_node[9,8] + star_edge[2,9] >= -2
 Partition_8_-_Star_edge_{2,_8}: - 2 partition[8,8]
   - circuit_partition_node[2,8] + star_edge[2,8] >= -2
 Partition_8_-_Star_edge_{2,_8}: - 2 partition[2,8]
   - circuit_partition_node[8,8] + star_edge[2,8] >= -2
 Partition_8_-_Star_edge_{2,_7}: - 2 partition[7,8]
   - circuit_partition_node[2,8] + star_edge[2,7] >= -2
 Partition_8_-_Star_edge_{2,_7}: - 2 partition[2,8]
   - circuit_partition_node[7,8] + star_edge[2,7] >= -2
 Partition_8_-_Star_edge_{2,_6}: - 2 partition[6,8]
   - circuit_partition_node[2,8] + star_edge[2,6] >= -2
 Partition_8_-_Star_edge_{2,_6}: - 2 partition[2,8]
   - circuit_partition_node[6,8] + star_edge[2,6] >= -2
 Partition_8_-_Star_edge_{2,_5}: - 2 partition[5,8]
   - circuit_partition_node[2,8] + star_edge[2,5] >= -2
 Partition_8_-_Star_edge_{2,_5}: - 2 partition[2,8]
   - circuit_partition_node[5,8] + star_edge[2,5] >= -2
 Partition_8_-_Star_edge_{2,_4}: - 2 partition[4,8]
   - circuit_partition_node[2,8] + star_edge[2,4] >= -2
 Partition_8_-_Star_edge_{2,_4}: - 2 partition[2,8]
   - circuit_partition_node[4,8] + star_edge[2,4] >= -2
 Partition_8_-_Star_edge_{2,_3}: - 2 partition[3,8]
   - circuit_partition_node[2,8] + star_edge[2,3] >= -2
 Partition_8_-_Star_edge_{2,_3}: - 2 partition[2,8]
   - circuit_partition_node[3,8] + star_edge[2,3] >= -2
 Partition_8_-_Star_edge_{6,_7}: - 2 partition[7,8]
   - circuit_partition_node[6,8] + star_edge[6,7] >= -2
 Partition_8_-_Star_edge_{6,_7}: - 2 partition[6,8]
   - circuit_partition_node[7,8] + star_edge[6,7] >= -2
 Partition_8_-_Star_edge_{6,_8}: - 2 partition[8,8]
   - circuit_partition_node[6,8] + star_edge[6,8] >= -2
 Partition_8_-_Star_edge_{6,_8}: - 2 partition[6,8]
   - circuit_partition_node[8,8] + star_edge[6,8] >= -2
 Partition_8_-_Star_edge_{6,_9}: - 2 partition[9,8]
   - circuit_partition_node[6,8] + star_edge[6,9] >= -2
 Partition_8_-_Star_edge_{6,_9}: - 2 partition[6,8]
   - circuit_partition_node[9,8] + star_edge[6,9] >= -2
 Partition_8_-_Star_edge_{1,_9}: - 2 partition[9,8]
   - circuit_partition_node[1,8] + star_edge[1,9] >= -2
 Partition_8_-_Star_edge_{1,_9}: - 2 partition[1,8]
   - circuit_partition_node[9,8] + star_edge[1,9] >= -2
 Partition_8_-_Star_edge_{1,_8}: - 2 partition[8,8]
   - circuit_partition_node[1,8] + star_edge[1,8] >= -2
 Partition_8_-_Star_edge_{1,_8}: - 2 partition[1,8]
   - circuit_partition_node[8,8] + star_edge[1,8] >= -2
 Partition_8_-_Star_edge_{1,_7}: - 2 partition[7,8]
   - circuit_partition_node[1,8] + star_edge[1,7] >= -2
 Partition_8_-_Star_edge_{1,_7}: - 2 partition[1,8]
   - circuit_partition_node[7,8] + star_edge[1,7] >= -2
 Partition_8_-_Star_edge_{1,_6}: - 2 partition[6,8]
   - circuit_partition_node[1,8] + star_edge[1,6] >= -2
 Partition_8_-_Star_edge_{1,_6}: - 2 partition[1,8]
   - circuit_partition_node[6,8] + star_edge[1,6] >= -2
 Partition_8_-_Star_edge_{1,_5}: - 2 partition[5,8]
   - circuit_partition_node[1,8] + star_edge[1,5] >= -2
 Partition_8_-_Star_edge_{1,_5}: - 2 partition[1,8]
   - circuit_partition_node[5,8] + star_edge[1,5] >= -2
 Partition_8_-_Star_edge_{1,_4}: - 2 partition[4,8]
   - circuit_partition_node[1,8] + star_edge[1,4] >= -2
 Partition_8_-_Star_edge_{1,_4}: - 2 partition[1,8]
   - circuit_partition_node[4,8] + star_edge[1,4] >= -2
 Partition_8_-_Star_edge_{1,_3}: - 2 partition[3,8]
   - circuit_partition_node[1,8] + star_edge[1,3] >= -2
 Partition_8_-_Star_edge_{1,_3}: - 2 partition[1,8]
   - circuit_partition_node[3,8] + star_edge[1,3] >= -2
 Partition_8_-_Star_edge_{1,_2}: - 2 partition[2,8]
   - circuit_partition_node[1,8] + star_edge[1,2] >= -2
 Partition_8_-_Star_edge_{1,_2}: - 2 partition[1,8]
   - circuit_partition_node[2,8] + star_edge[1,2] >= -2
 Partition_8_-_Star_edge_{7,_8}: - 2 partition[8,8]
   - circuit_partition_node[7,8] + star_edge[7,8] >= -2
 Partition_8_-_Star_edge_{7,_8}: - 2 partition[7,8]
   - circuit_partition_node[8,8] + star_edge[7,8] >= -2
 Partition_8_-_Star_edge_{7,_9}: - 2 partition[9,8]
   - circuit_partition_node[7,8] + star_edge[7,9] >= -2
 Partition_8_-_Star_edge_{7,_9}: - 2 partition[7,8]
   - circuit_partition_node[9,8] + star_edge[7,9] >= -2
 Partition_8_-_Star_edge_{0,_9}: - 2 partition[9,8]
   - circuit_partition_node[0,8] + star_edge[0,9] >= -2
 Partition_8_-_Star_edge_{0,_9}: - 2 partition[0,8]
   - circuit_partition_node[9,8] + star_edge[0,9] >= -2
 Partition_8_-_Star_edge_{0,_8}: - 2 partition[8,8]
   - circuit_partition_node[0,8] + star_edge[0,8] >= -2
 Partition_8_-_Star_edge_{0,_8}: - 2 partition[0,8]
   - circuit_partition_node[8,8] + star_edge[0,8] >= -2
 Partition_8_-_Star_edge_{0,_7}: - 2 partition[7,8]
   - circuit_partition_node[0,8] + star_edge[0,7] >= -2
 Partition_8_-_Star_edge_{0,_7}: - 2 partition[0,8]
   - circuit_partition_node[7,8] + star_edge[0,7] >= -2
 Partition_8_-_Star_edge_{0,_6}: - 2 partition[6,8]
   - circuit_partition_node[0,8] + star_edge[0,6] >= -2
 Partition_8_-_Star_edge_{0,_6}: - 2 partition[0,8]
   - circuit_partition_node[6,8] + star_edge[0,6] >= -2
 Partition_8_-_Star_edge_{0,_5}: - 2 partition[5,8]
   - circuit_partition_node[0,8] + star_edge[0,5] >= -2
 Partition_8_-_Star_edge_{0,_5}: - 2 partition[0,8]
   - circuit_partition_node[5,8] + star_edge[0,5] >= -2
 Partition_8_-_Star_edge_{0,_4}: - 2 partition[4,8]
   - circuit_partition_node[0,8] + star_edge[0,4] >= -2
 Partition_8_-_Star_edge_{0,_4}: - 2 partition[0,8]
   - circuit_partition_node[4,8] + star_edge[0,4] >= -2
 Partition_8_-_Star_edge_{0,_3}: - 2 partition[3,8]
   - circuit_partition_node[0,8] + star_edge[0,3] >= -2
 Partition_8_-_Star_edge_{0,_3}: - 2 partition[0,8]
   - circuit_partition_node[3,8] + star_edge[0,3] >= -2
 Partition_8_-_Star_edge_{0,_2}: - 2 partition[2,8]
   - circuit_partition_node[0,8] + star_edge[0,2] >= -2
 Partition_8_-_Star_edge_{0,_2}: - 2 partition[0,8]
   - circuit_partition_node[2,8] + star_edge[0,2] >= -2
 Partition_8_-_Star_edge_{0,_1}: - 2 partition[1,8]
   - circuit_partition_node[0,8] + star_edge[0,1] >= -2
 Partition_8_-_Star_edge_{0,_1}: - 2 partition[0,8]
   - circuit_partition_node[1,8] + star_edge[0,1] >= -2
 Partition_8_-_Star_edge_{8,_9}: - 2 partition[9,8]
   - circuit_partition_node[8,8] + star_edge[8,9] >= -2
 Partition_8_-_Star_edge_{8,_9}: - 2 partition[8,8]
   - circuit_partition_node[9,8] + star_edge[8,9] >= -2
 Partition_9_-_Star_edge_{4,_5}: - 2 partition[5,9]
   - circuit_partition_node[4,9] + star_edge[4,5] >= -2
 Partition_9_-_Star_edge_{4,_5}: - 2 partition[4,9]
   - circuit_partition_node[5,9] + star_edge[4,5] >= -2
 Partition_9_-_Star_edge_{4,_6}: - 2 partition[6,9]
   - circuit_partition_node[4,9] + star_edge[4,6] >= -2
 Partition_9_-_Star_edge_{4,_6}: - 2 partition[4,9]
   - circuit_partition_node[6,9] + star_edge[4,6] >= -2
 Partition_9_-_Star_edge_{4,_7}: - 2 partition[7,9]
   - circuit_partition_node[4,9] + star_edge[4,7] >= -2
 Partition_9_-_Star_edge_{4,_7}: - 2 partition[4,9]
   - circuit_partition_node[7,9] + star_edge[4,7] >= -2
 Partition_9_-_Star_edge_{4,_8}: - 2 partition[8,9]
   - circuit_partition_node[4,9] + star_edge[4,8] >= -2
 Partition_9_-_Star_edge_{4,_8}: - 2 partition[4,9]
   - circuit_partition_node[8,9] + star_edge[4,8] >= -2
 Partition_9_-_Star_edge_{4,_9}: - 2 partition[9,9]
   - circuit_partition_node[4,9] + star_edge[4,9] >= -2
 Partition_9_-_Star_edge_{4,_9}: - 2 partition[4,9]
   - circuit_partition_node[9,9] + star_edge[4,9] >= -2
 Partition_9_-_Star_edge_{3,_9}: - 2 partition[9,9]
   - circuit_partition_node[3,9] + star_edge[3,9] >= -2
 Partition_9_-_Star_edge_{3,_9}: - 2 partition[3,9]
   - circuit_partition_node[9,9] + star_edge[3,9] >= -2
 Partition_9_-_Star_edge_{3,_8}: - 2 partition[8,9]
   - circuit_partition_node[3,9] + star_edge[3,8] >= -2
 Partition_9_-_Star_edge_{3,_8}: - 2 partition[3,9]
   - circuit_partition_node[8,9] + star_edge[3,8] >= -2
 Partition_9_-_Star_edge_{3,_7}: - 2 partition[7,9]
   - circuit_partition_node[3,9] + star_edge[3,7] >= -2
 Partition_9_-_Star_edge_{3,_7}: - 2 partition[3,9]
   - circuit_partition_node[7,9] + star_edge[3,7] >= -2
 Partition_9_-_Star_edge_{3,_6}: - 2 partition[6,9]
   - circuit_partition_node[3,9] + star_edge[3,6] >= -2
 Partition_9_-_Star_edge_{3,_6}: - 2 partition[3,9]
   - circuit_partition_node[6,9] + star_edge[3,6] >= -2
 Partition_9_-_Star_edge_{3,_5}: - 2 partition[5,9]
   - circuit_partition_node[3,9] + star_edge[3,5] >= -2
 Partition_9_-_Star_edge_{3,_5}: - 2 partition[3,9]
   - circuit_partition_node[5,9] + star_edge[3,5] >= -2
 Partition_9_-_Star_edge_{3,_4}: - 2 partition[4,9]
   - circuit_partition_node[3,9] + star_edge[3,4] >= -2
 Partition_9_-_Star_edge_{3,_4}: - 2 partition[3,9]
   - circuit_partition_node[4,9] + star_edge[3,4] >= -2
 Partition_9_-_Star_edge_{5,_6}: - 2 partition[6,9]
   - circuit_partition_node[5,9] + star_edge[5,6] >= -2
 Partition_9_-_Star_edge_{5,_6}: - 2 partition[5,9]
   - circuit_partition_node[6,9] + star_edge[5,6] >= -2
 Partition_9_-_Star_edge_{5,_7}: - 2 partition[7,9]
   - circuit_partition_node[5,9] + star_edge[5,7] >= -2
 Partition_9_-_Star_edge_{5,_7}: - 2 partition[5,9]
   - circuit_partition_node[7,9] + star_edge[5,7] >= -2
 Partition_9_-_Star_edge_{5,_8}: - 2 partition[8,9]
   - circuit_partition_node[5,9] + star_edge[5,8] >= -2
 Partition_9_-_Star_edge_{5,_8}: - 2 partition[5,9]
   - circuit_partition_node[8,9] + star_edge[5,8] >= -2
 Partition_9_-_Star_edge_{5,_9}: - 2 partition[9,9]
   - circuit_partition_node[5,9] + star_edge[5,9] >= -2
 Partition_9_-_Star_edge_{5,_9}: - 2 partition[5,9]
   - circuit_partition_node[9,9] + star_edge[5,9] >= -2
 Partition_9_-_Star_edge_{2,_9}: - 2 partition[9,9]
   - circuit_partition_node[2,9] + star_edge[2,9] >= -2
 Partition_9_-_Star_edge_{2,_9}: - 2 partition[2,9]
   - circuit_partition_node[9,9] + star_edge[2,9] >= -2
 Partition_9_-_Star_edge_{2,_8}: - 2 partition[8,9]
   - circuit_partition_node[2,9] + star_edge[2,8] >= -2
 Partition_9_-_Star_edge_{2,_8}: - 2 partition[2,9]
   - circuit_partition_node[8,9] + star_edge[2,8] >= -2
 Partition_9_-_Star_edge_{2,_7}: - 2 partition[7,9]
   - circuit_partition_node[2,9] + star_edge[2,7] >= -2
 Partition_9_-_Star_edge_{2,_7}: - 2 partition[2,9]
   - circuit_partition_node[7,9] + star_edge[2,7] >= -2
 Partition_9_-_Star_edge_{2,_6}: - 2 partition[6,9]
   - circuit_partition_node[2,9] + star_edge[2,6] >= -2
 Partition_9_-_Star_edge_{2,_6}: - 2 partition[2,9]
   - circuit_partition_node[6,9] + star_edge[2,6] >= -2
 Partition_9_-_Star_edge_{2,_5}: - 2 partition[5,9]
   - circuit_partition_node[2,9] + star_edge[2,5] >= -2
 Partition_9_-_Star_edge_{2,_5}: - 2 partition[2,9]
   - circuit_partition_node[5,9] + star_edge[2,5] >= -2
 Partition_9_-_Star_edge_{2,_4}: - 2 partition[4,9]
   - circuit_partition_node[2,9] + star_edge[2,4] >= -2
 Partition_9_-_Star_edge_{2,_4}: - 2 partition[2,9]
   - circuit_partition_node[4,9] + star_edge[2,4] >= -2
 Partition_9_-_Star_edge_{2,_3}: - 2 partition[3,9]
   - circuit_partition_node[2,9] + star_edge[2,3] >= -2
 Partition_9_-_Star_edge_{2,_3}: - 2 partition[2,9]
   - circuit_partition_node[3,9] + star_edge[2,3] >= -2
 Partition_9_-_Star_edge_{6,_7}: - 2 partition[7,9]
   - circuit_partition_node[6,9] + star_edge[6,7] >= -2
 Partition_9_-_Star_edge_{6,_7}: - 2 partition[6,9]
   - circuit_partition_node[7,9] + star_edge[6,7] >= -2
 Partition_9_-_Star_edge_{6,_8}: - 2 partition[8,9]
   - circuit_partition_node[6,9] + star_edge[6,8] >= -2
 Partition_9_-_Star_edge_{6,_8}: - 2 partition[6,9]
   - circuit_partition_node[8,9] + star_edge[6,8] >= -2
 Partition_9_-_Star_edge_{6,_9}: - 2 partition[9,9]
   - circuit_partition_node[6,9] + star_edge[6,9] >= -2
 Partition_9_-_Star_edge_{6,_9}: - 2 partition[6,9]
   - circuit_partition_node[9,9] + star_edge[6,9] >= -2
 Partition_9_-_Star_edge_{1,_9}: - 2 partition[9,9]
   - circuit_partition_node[1,9] + star_edge[1,9] >= -2
 Partition_9_-_Star_edge_{1,_9}: - 2 partition[1,9]
   - circuit_partition_node[9,9] + star_edge[1,9] >= -2
 Partition_9_-_Star_edge_{1,_8}: - 2 partition[8,9]
   - circuit_partition_node[1,9] + star_edge[1,8] >= -2
 Partition_9_-_Star_edge_{1,_8}: - 2 partition[1,9]
   - circuit_partition_node[8,9] + star_edge[1,8] >= -2
 Partition_9_-_Star_edge_{1,_7}: - 2 partition[7,9]
   - circuit_partition_node[1,9] + star_edge[1,7] >= -2
 Partition_9_-_Star_edge_{1,_7}: - 2 partition[1,9]
   - circuit_partition_node[7,9] + star_edge[1,7] >= -2
 Partition_9_-_Star_edge_{1,_6}: - 2 partition[6,9]
   - circuit_partition_node[1,9] + star_edge[1,6] >= -2
 Partition_9_-_Star_edge_{1,_6}: - 2 partition[1,9]
   - circuit_partition_node[6,9] + star_edge[1,6] >= -2
 Partition_9_-_Star_edge_{1,_5}: - 2 partition[5,9]
   - circuit_partition_node[1,9] + star_edge[1,5] >= -2
 Partition_9_-_Star_edge_{1,_5}: - 2 partition[1,9]
   - circuit_partition_node[5,9] + star_edge[1,5] >= -2
 Partition_9_-_Star_edge_{1,_4}: - 2 partition[4,9]
   - circuit_partition_node[1,9] + star_edge[1,4] >= -2
 Partition_9_-_Star_edge_{1,_4}: - 2 partition[1,9]
   - circuit_partition_node[4,9] + star_edge[1,4] >= -2
 Partition_9_-_Star_edge_{1,_3}: - 2 partition[3,9]
   - circuit_partition_node[1,9] + star_edge[1,3] >= -2
 Partition_9_-_Star_edge_{1,_3}: - 2 partition[1,9]
   - circuit_partition_node[3,9] + star_edge[1,3] >= -2
 Partition_9_-_Star_edge_{1,_2}: - 2 partition[2,9]
   - circuit_partition_node[1,9] + star_edge[1,2] >= -2
 Partition_9_-_Star_edge_{1,_2}: - 2 partition[1,9]
   - circuit_partition_node[2,9] + star_edge[1,2] >= -2
 Partition_9_-_Star_edge_{7,_8}: - 2 partition[8,9]
   - circuit_partition_node[7,9] + star_edge[7,8] >= -2
 Partition_9_-_Star_edge_{7,_8}: - 2 partition[7,9]
   - circuit_partition_node[8,9] + star_edge[7,8] >= -2
 Partition_9_-_Star_edge_{7,_9}: - 2 partition[9,9]
   - circuit_partition_node[7,9] + star_edge[7,9] >= -2
 Partition_9_-_Star_edge_{7,_9}: - 2 partition[7,9]
   - circuit_partition_node[9,9] + star_edge[7,9] >= -2
 Partition_9_-_Star_edge_{0,_9}: - 2 partition[9,9]
   - circuit_partition_node[0,9] + star_edge[0,9] >= -2
 Partition_9_-_Star_edge_{0,_9}: - 2 partition[0,9]
   - circuit_partition_node[9,9] + star_edge[0,9] >= -2
 Partition_9_-_Star_edge_{0,_8}: - 2 partition[8,9]
   - circuit_partition_node[0,9] + star_edge[0,8] >= -2
 Partition_9_-_Star_edge_{0,_8}: - 2 partition[0,9]
   - circuit_partition_node[8,9] + star_edge[0,8] >= -2
 Partition_9_-_Star_edge_{0,_7}: - 2 partition[7,9]
   - circuit_partition_node[0,9] + star_edge[0,7] >= -2
 Partition_9_-_Star_edge_{0,_7}: - 2 partition[0,9]
   - circuit_partition_node[7,9] + star_edge[0,7] >= -2
 Partition_9_-_Star_edge_{0,_6}: - 2 partition[6,9]
   - circuit_partition_node[0,9] + star_edge[0,6] >= -2
 Partition_9_-_Star_edge_{0,_6}: - 2 partition[0,9]
   - circuit_partition_node[6,9] + star_edge[0,6] >= -2
 Partition_9_-_Star_edge_{0,_5}: - 2 partition[5,9]
   - circuit_partition_node[0,9] + star_edge[0,5] >= -2
 Partition_9_-_Star_edge_{0,_5}: - 2 partition[0,9]
   - circuit_partition_node[5,9] + star_edge[0,5] >= -2
 Partition_9_-_Star_edge_{0,_4}: - 2 partition[4,9]
   - circuit_partition_node[0,9] + star_edge[0,4] >= -2
 Partition_9_-_Star_edge_{0,_4}: - 2 partition[0,9]
   - circuit_partition_node[4,9] + star_edge[0,4] >= -2
 Partition_9_-_Star_edge_{0,_3}: - 2 partition[3,9]
   - circuit_partition_node[0,9] + star_edge[0,3] >= -2
 Partition_9_-_Star_edge_{0,_3}: - 2 partition[0,9]
   - circuit_partition_node[3,9] + star_edge[0,3] >= -2
 Partition_9_-_Star_edge_{0,_2}: - 2 partition[2,9]
   - circuit_partition_node[0,9] + star_edge[0,2] >= -2
 Partition_9_-_Star_edge_{0,_2}: - 2 partition[0,9]
   - circuit_partition_node[2,9] + star_edge[0,2] >= -2
 Partition_9_-_Star_edge_{0,_1}: - 2 partition[1,9]
   - circuit_partition_node[0,9] + star_edge[0,1] >= -2
 Partition_9_-_Star_edge_{0,_1}: - 2 partition[0,9]
   - circuit_partition_node[1,9] + star_edge[0,1] >= -2
 Partition_9_-_Star_edge_{8,_9}: - 2 partition[9,9]
   - circuit_partition_node[8,9] + star_edge[8,9] >= -2
 Partition_9_-_Star_edge_{8,_9}: - 2 partition[8,9]
   - circuit_partition_node[9,9] + star_edge[8,9] >= -2
 Partition_0_-_Exactly_one_circuit_node: - partition_used[0]
   + circuit_partition_node[9,0] + circuit_partition_node[8,0]
   + circuit_partition_node[7,0] + circuit_partition_node[6,0]
   + circuit_partition_node[5,0] + circuit_partition_node[4,0]
   + circuit_partition_node[3,0] + circuit_partition_node[2,0]
   + circuit_partition_node[1,0] + circuit_partition_node[0,0] = 0
 Partition_0_-_Circuit_node_9_validity: - partition[9,0]
   + circuit_partition_node[9,0] <= 0
 Partition_0_-_Circuit_node_8_validity: - partition[8,0]
   + circuit_partition_node[8,0] <= 0
 Partition_0_-_Circuit_node_7_validity: - partition[7,0]
   + circuit_partition_node[7,0] <= 0
 Partition_0_-_Circuit_node_6_validity: - partition[6,0]
   + circuit_partition_node[6,0] <= 0
 Partition_0_-_Circuit_node_5_validity: - partition[5,0]
   + circuit_partition_node[5,0] <= 0
 Partition_0_-_Circuit_node_4_validity: - partition[4,0]
   + circuit_partition_node[4,0] <= 0
 Partition_0_-_Circuit_node_3_validity: - partition[3,0]
   + circuit_partition_node[3,0] <= 0
 Partition_0_-_Circuit_node_2_validity: - partition[2,0]
   + circuit_partition_node[2,0] <= 0
 Partition_0_-_Circuit_node_1_validity: - partition[1,0]
   + circuit_partition_node[1,0] <= 0
 Partition_0_-_Circuit_node_0_validity: - partition[0,0]
   + circuit_partition_node[0,0] <= 0
 Partition_1_-_Exactly_one_circuit_node: - partition_used[1]
   + circuit_partition_node[9,1] + circuit_partition_node[8,1]
   + circuit_partition_node[7,1] + circuit_partition_node[6,1]
   + circuit_partition_node[5,1] + circuit_partition_node[4,1]
   + circuit_partition_node[3,1] + circuit_partition_node[2,1]
   + circuit_partition_node[1,1] + circuit_partition_node[0,1] = 0
 Partition_1_-_Circuit_node_9_validity: - partition[9,1]
   + circuit_partition_node[9,1] <= 0
 Partition_1_-_Circuit_node_8_validity: - partition[8,1]
   + circuit_partition_node[8,1] <= 0
 Partition_1_-_Circuit_node_7_validity: - partition[7,1]
   + circuit_partition_node[7,1] <= 0
 Partition_1_-_Circuit_node_6_validity: - partition[6,1]
   + circuit_partition_node[6,1] <= 0
 Partition_1_-_Circuit_node_5_validity: - partition[5,1]
   + circuit_partition_node[5,1] <= 0
 Partition_1_-_Circuit_node_4_validity: - partition[4,1]
   + circuit_partition_node[4,1] <= 0
 Partition_1_-_Circuit_node_3_validity: - partition[3,1]
   + circuit_partition_node[3,1] <= 0
 Partition_1_-_Circuit_node_2_validity: - partition[2,1]
   + circuit_partition_node[2,1] <= 0
 Partition_1_-_Circuit_node_1_validity: - partition[1,1]
   + circuit_partition_node[1,1] <= 0
 Partition_1_-_Circuit_node_0_validity: - partition[0,1]
   + circuit_partition_node[0,1] <= 0
 Partition_2_-_Exactly_one_circuit_node: - partition_used[2]
   + circuit_partition_node[9,2] + circuit_partition_node[8,2]
   + circuit_partition_node[7,2] + circuit_partition_node[6,2]
   + circuit_partition_node[5,2] + circuit_partition_node[4,2]
   + circuit_partition_node[3,2] + circuit_partition_node[2,2]
   + circuit_partition_node[1,2] + circuit_partition_node[0,2] = 0
 Partition_2_-_Circuit_node_9_validity: - partition[9,2]
   + circuit_partition_node[9,2] <= 0
 Partition_2_-_Circuit_node_8_validity: - partition[8,2]
   + circuit_partition_node[8,2] <= 0
 Partition_2_-_Circuit_node_7_validity: - partition[7,2]
   + circuit_partition_node[7,2] <= 0
 Partition_2_-_Circuit_node_6_validity: - partition[6,2]
   + circuit_partition_node[6,2] <= 0
 Partition_2_-_Circuit_node_5_validity: - partition[5,2]
   + circuit_partition_node[5,2] <= 0
 Partition_2_-_Circuit_node_4_validity: - partition[4,2]
   + circuit_partition_node[4,2] <= 0
 Partition_2_-_Circuit_node_3_validity: - partition[3,2]
   + circuit_partition_node[3,2] <= 0
 Partition_2_-_Circuit_node_2_validity: - partition[2,2]
   + circuit_partition_node[2,2] <= 0
 Partition_2_-_Circuit_node_1_validity: - partition[1,2]
   + circuit_partition_node[1,2] <= 0
 Partition_2_-_Circuit_node_0_validity: - partition[0,2]
   + circuit_partition_node[0,2] <= 0
 Partition_3_-_Exactly_one_circuit_node: - partition_used[3]
   + circuit_partition_node[9,3] + circuit_partition_node[8,3]
   + circuit_partition_node[7,3] + circuit_partition_node[6,3]
   + circuit_partition_node[5,3] + circuit_partition_node[4,3]
   + circuit_partition_node[3,3] + circuit_partition_node[2,3]
   + circuit_partition_node[1,3] + circuit_partition_node[0,3] = 0
 Partition_3_-_Circuit_node_9_validity: - partition[9,3]
   + circuit_partition_node[9,3] <= 0
 Partition_3_-_Circuit_node_8_validity: - partition[8,3]
   + circuit_partition_node[8,3] <= 0
 Partition_3_-_Circuit_node_7_validity: - partition[7,3]
   + circuit_partition_node[7,3] <= 0
 Partition_3_-_Circuit_node_6_validity: - partition[6,3]
   + circuit_partition_node[6,3] <= 0
 Partition_3_-_Circuit_node_5_validity: - partition[5,3]
   + circuit_partition_node[5,3] <= 0
 Partition_3_-_Circuit_node_4_validity: - partition[4,3]
   + circuit_partition_node[4,3] <= 0
 Partition_3_-_Circuit_node_3_validity: - partition[3,3]
   + circuit_partition_node[3,3] <= 0
 Partition_3_-_Circuit_node_2_validity: - partition[2,3]
   + circuit_partition_node[2,3] <= 0
 Partition_3_-_Circuit_node_1_validity: - partition[1,3]
   + circuit_partition_node[1,3] <= 0
 Partition_3_-_Circuit_node_0_validity: - partition[0,3]
   + circuit_partition_node[0,3] <= 0
 Partition_4_-_Exactly_one_circuit_node: - partition_used[4]
   + circuit_partition_node[9,4] + circuit_partition_node[8,4]
   + circuit_partition_node[7,4] + circuit_partition_node[6,4]
   + circuit_partition_node[5,4] + circuit_partition_node[4,4]
   + circuit_partition_node[3,4] + circuit_partition_node[2,4]
   + circuit_partition_node[1,4] + circuit_partition_node[0,4] = 0
 Partition_4_-_Circuit_node_9_validity: - partition[9,4]
   + circuit_partition_node[9,4] <= 0
 Partition_4_-_Circuit_node_8_validity: - partition[8,4]
   + circuit_partition_node[8,4] <= 0
 Partition_4_-_Circuit_node_7_validity: - partition[7,4]
   + circuit_partition_node[7,4] <= 0
 Partition_4_-_Circuit_node_6_validity: - partition[6,4]
   + circuit_partition_node[6,4] <= 0
 Partition_4_-_Circuit_node_5_validity: - partition[5,4]
   + circuit_partition_node[5,4] <= 0
 Partition_4_-_Circuit_node_4_validity: - partition[4,4]
   + circuit_partition_node[4,4] <= 0
 Partition_4_-_Circuit_node_3_validity: - partition[3,4]
   + circuit_partition_node[3,4] <= 0
 Partition_4_-_Circuit_node_2_validity: - partition[2,4]
   + circuit_partition_node[2,4] <= 0
 Partition_4_-_Circuit_node_1_validity: - partition[1,4]
   + circuit_partition_node[1,4] <= 0
 Partition_4_-_Circuit_node_0_validity: - partition[0,4]
   + circuit_partition_node[0,4] <= 0
 Partition_5_-_Exactly_one_circuit_node: - partition_used[5]
   + circuit_partition_node[9,5] + circuit_partition_node[8,5]
   + circuit_partition_node[7,5] + circuit_partition_node[6,5]
   + circuit_partition_node[5,5] + circuit_partition_node[4,5]
   + circuit_partition_node[3,5] + circuit_partition_node[2,5]
   + circuit_partition_node[1,5] + circuit_partition_node[0,5] = 0
 Partition_5_-_Circuit_node_9_validity: - partition[9,5]
   + circuit_partition_node[9,5] <= 0
 Partition_5_-_Circuit_node_8_validity: - partition[8,5]
   + circuit_partition_node[8,5] <= 0
 Partition_5_-_Circuit_node_7_validity: - partition[7,5]
   + circuit_partition_node[7,5] <= 0
 Partition_5_-_Circuit_node_6_validity: - partition[6,5]
   + circuit_partition_node[6,5] <= 0
 Partition_5_-_Circuit_node_5_validity: - partition[5,5]
   + circuit_partition_node[5,5] <= 0
 Partition_5_-_Circuit_node_4_validity: - partition[4,5]
   + circuit_partition_node[4,5] <= 0
 Partition_5_-_Circuit_node_3_validity: - partition[3,5]
   + circuit_partition_node[3,5] <= 0
 Partition_5_-_Circuit_node_2_validity: - partition[2,5]
   + circuit_partition_node[2,5] <= 0
 Partition_5_-_Circuit_node_1_validity: - partition[1,5]
   + circuit_partition_node[1,5] <= 0
 Partition_5_-_Circuit_node_0_validity: - partition[0,5]
   + circuit_partition_node[0,5] <= 0
 Partition_6_-_Exactly_one_circuit_node: - partition_used[6]
   + circuit_partition_node[9,6] + circuit_partition_node[8,6]
   + circuit_partition_node[7,6] + circuit_partition_node[6,6]
   + circuit_partition_node[5,6] + circuit_partition_node[4,6]
   + circuit_partition_node[3,6] + circuit_partition_node[2,6]
   + circuit_partition_node[1,6] + circuit_partition_node[0,6] = 0
 Partition_6_-_Circuit_node_9_validity: - partition[9,6]
   + circuit_partition_node[9,6] <= 0
 Partition_6_-_Circuit_node_8_validity: - partition[8,6]
   + circuit_partition_node[8,6] <= 0
 Partition_6_-_Circuit_node_7_validity: - partition[7,6]
   + circuit_partition_node[7,6] <= 0
 Partition_6_-_Circuit_node_6_validity: - partition[6,6]
   + circuit_partition_node[6,6] <= 0
 Partition_6_-_Circuit_node_5_validity: - partition[5,6]
   + circuit_partition_node[5,6] <= 0
 Partition_6_-_Circuit_node_4_validity: - partition[4,6]
   + circuit_partition_node[4,6] <= 0
 Partition_6_-_Circuit_node_3_validity: - partition[3,6]
   + circuit_partition_node[3,6] <= 0
 Partition_6_-_Circuit_node_2_validity: - partition[2,6]
   + circuit_partition_node[2,6] <= 0
 Partition_6_-_Circuit_node_1_validity: - partition[1,6]
   + circuit_partition_node[1,6] <= 0
 Partition_6_-_Circuit_node_0_validity: - partition[0,6]
   + circuit_partition_node[0,6] <= 0
 Partition_7_-_Exactly_one_circuit_node: - partition_used[7]
   + circuit_partition_node[9,7] + circuit_partition_node[8,7]
   + circuit_partition_node[7,7] + circuit_partition_node[6,7]
   + circuit_partition_node[5,7] + circuit_partition_node[4,7]
   + circuit_partition_node[3,7] + circuit_partition_node[2,7]
   + circuit_partition_node[1,7] + circuit_partition_node[0,7] = 0
 Partition_7_-_Circuit_node_9_validity: - partition[9,7]
   + circuit_partition_node[9,7] <= 0
 Partition_7_-_Circuit_node_8_validity: - partition[8,7]
   + circuit_partition_node[8,7] <= 0
 Partition_7_-_Circuit_node_7_validity: - partition[7,7]
   + circuit_partition_node[7,7] <= 0
 Partition_7_-_Circuit_node_6_validity: - partition[6,7]
   + circuit_partition_node[6,7] <= 0
 Partition_7_-_Circuit_node_5_validity: - partition[5,7]
   + circuit_partition_node[5,7] <= 0
 Partition_7_-_Circuit_node_4_validity: - partition[4,7]
   + circuit_partition_node[4,7] <= 0
 Partition_7_-_Circuit_node_3_validity: - partition[3,7]
   + circuit_partition_node[3,7] <= 0
 Partition_7_-_Circuit_node_2_validity: - partition[2,7]
   + circuit_partition_node[2,7] <= 0
 Partition_7_-_Circuit_node_1_validity: - partition[1,7]
   + circuit_partition_node[1,7] <= 0
 Partition_7_-_Circuit_node_0_validity: - partition[0,7]
   + circuit_partition_node[0,7] <= 0
 Partition_8_-_Exactly_one_circuit_node: - partition_used[8]
   + circuit_partition_node[9,8] + circuit_partition_node[8,8]
   + circuit_partition_node[7,8] + circuit_partition_node[6,8]
   + circuit_partition_node[5,8] + circuit_partition_node[4,8]
   + circuit_partition_node[3,8] + circuit_partition_node[2,8]
   + circuit_partition_node[1,8] + circuit_partition_node[0,8] = 0
 Partition_8_-_Circuit_node_9_validity: - partition[9,8]
   + circuit_partition_node[9,8] <= 0
 Partition_8_-_Circuit_node_8_validity: - partition[8,8]
   + circuit_partition_node[8,8] <= 0
 Partition_8_-_Circuit_node_7_validity: - partition[7,8]
   + circuit_partition_node[7,8] <= 0
 Partition_8_-_Circuit_node_6_validity: - partition[6,8]
   + circuit_partition_node[6,8] <= 0
 Partition_8_-_Circuit_node_5_validity: - partition[5,8]
   + circuit_partition_node[5,8] <= 0
 Partition_8_-_Circuit_node_4_validity: - partition[4,8]
   + circuit_partition_node[4,8] <= 0
 Partition_8_-_Circuit_node_3_validity: - partition[3,8]
   + circuit_partition_node[3,8] <= 0
 Partition_8_-_Circuit_node_2_validity: - partition[2,8]
   + circuit_partition_node[2,8] <= 0
 Partition_8_-_Circuit_node_1_validity: - partition[1,8]
   + circuit_partition_node[1,8] <= 0
 Partition_8_-_Circuit_node_0_validity: - partition[0,8]
   + circuit_partition_node[0,8] <= 0
 Partition_9_-_Exactly_one_circuit_node: - partition_used[9]
   + circuit_partition_node[9,9] + circuit_partition_node[8,9]
   + circuit_partition_node[7,9] + circuit_partition_node[6,9]
   + circuit_partition_node[5,9] + circuit_partition_node[4,9]
   + circuit_partition_node[3,9] + circuit_partition_node[2,9]
   + circuit_partition_node[1,9] + circuit_partition_node[0,9] = 0
 Partition_9_-_Circuit_node_9_validity: - partition[9,9]
   + circuit_partition_node[9,9] <= 0
 Partition_9_-_Circuit_node_8_validity: - partition[8,9]
   + circuit_partition_node[8,9] <= 0
 Partition_9_-_Circuit_node_7_validity: - partition[7,9]
   + circuit_partition_node[7,9] <= 0
 Partition_9_-_Circuit_node_6_validity: - partition[6,9]
   + circuit_partition_node[6,9] <= 0
 Partition_9_-_Circuit_node_5_validity: - partition[5,9]
   + circuit_partition_node[5,9] <= 0
 Partition_9_-_Circuit_node_4_validity: - partition[4,9]
   + circuit_partition_node[4,9] <= 0
 Partition_9_-_Circuit_node_3_validity: - partition[3,9]
   + circuit_partition_node[3,9] <= 0
 Partition_9_-_Circuit_node_2_validity: - partition[2,9]
   + circuit_partition_node[2,9] <= 0
 Partition_9_-_Circuit_node_1_validity: - partition[1,9]
   + circuit_partition_node[1,9] <= 0
 Partition_9_-_Circuit_node_0_validity: - partition[0,9]
   + circuit_partition_node[0,9] <= 0
 R1140: - circuit_partition_node[9,0] + circuit_node[9] >= 0
 R1141: - circuit_partition_node[9,1] + circuit_node[9] >= 0
 R1142: - circuit_partition_node[9,2] + circuit_node[9] >= 0
 R1143: - circuit_partition_node[9,3] + circuit_node[9] >= 0
 R1144: - circuit_partition_node[9,4] + circuit_node[9] >= 0
 R1145: - circuit_partition_node[9,5] + circuit_node[9] >= 0
 R1146: - circuit_partition_node[9,6] + circuit_node[9] >= 0
 R1147: - circuit_partition_node[9,7] + circuit_node[9] >= 0
 R1148: - circuit_partition_node[9,8] + circuit_node[9] >= 0
 R1149: - circuit_partition_node[9,9] + circuit_node[9] >= 0
 R1150: - circuit_partition_node[8,0] + circuit_node[8] >= 0
 R1151: - circuit_partition_node[8,1] + circuit_node[8] >= 0
 R1152: - circuit_partition_node[8,2] + circuit_node[8] >= 0
 R1153: - circuit_partition_node[8,3] + circuit_node[8] >= 0
 R1154: - circuit_partition_node[8,4] + circuit_node[8] >= 0
 R1155: - circuit_partition_node[8,5] + circuit_node[8] >= 0
 R1156: - circuit_partition_node[8,6] + circuit_node[8] >= 0
 R1157: - circuit_partition_node[8,7] + circuit_node[8] >= 0
 R1158: - circuit_partition_node[8,8] + circuit_node[8] >= 0
 R1159: - circuit_partition_node[8,9] + circuit_node[8] >= 0
 R1160: - circuit_partition_node[7,0] + circuit_node[7] >= 0
 R1161: - circuit_partition_node[7,1] + circuit_node[7] >= 0
 R1162: - circuit_partition_node[7,2] + circuit_node[7] >= 0
 R1163: - circuit_partition_node[7,3] + circuit_node[7] >= 0
 R1164: - circuit_partition_node[7,4] + circuit_node[7] >= 0
 R1165: - circuit_partition_node[7,5] + circuit_node[7] >= 0
 R1166: - circuit_partition_node[7,6] + circuit_node[7] >= 0
 R1167: - circuit_partition_node[7,7] + circuit_node[7] >= 0
 R1168: - circuit_partition_node[7,8] + circuit_node[7] >= 0
 R1169: - circuit_partition_node[7,9] + circuit_node[7] >= 0
 R1170: - circuit_partition_node[6,0] + circuit_node[6] >= 0
 R1171: - circuit_partition_node[6,1] + circuit_node[6] >= 0
 R1172: - circuit_partition_node[6,2] + circuit_node[6] >= 0
 R1173: - circuit_partition_node[6,3] + circuit_node[6] >= 0
 R1174: - circuit_partition_node[6,4] + circuit_node[6] >= 0
 R1175: - circuit_partition_node[6,5] + circuit_node[6] >= 0
 R1176: - circuit_partition_node[6,6] + circuit_node[6] >= 0
 R1177: - circuit_partition_node[6,7] + circuit_node[6] >= 0
 R1178: - circuit_partition_node[6,8] + circuit_node[6] >= 0
 R1179: - circuit_partition_node[6,9] + circuit_node[6] >= 0
 R1180: - circuit_partition_node[5,0] + circuit_node[5] >= 0
 R1181: - circuit_partition_node[5,1] + circuit_node[5] >= 0
 R1182: - circuit_partition_node[5,2] + circuit_node[5] >= 0
 R1183: - circuit_partition_node[5,3] + circuit_node[5] >= 0
 R1184: - circuit_partition_node[5,4] + circuit_node[5] >= 0
 R1185: - circuit_partition_node[5,5] + circuit_node[5] >= 0
 R1186: - circuit_partition_node[5,6] + circuit_node[5] >= 0
 R1187: - circuit_partition_node[5,7] + circuit_node[5] >= 0
 R1188: - circuit_partition_node[5,8] + circuit_node[5] >= 0
 R1189: - circuit_partition_node[5,9] + circuit_node[5] >= 0
 R1190: - circuit_partition_node[4,0] + circuit_node[4] >= 0
 R1191: - circuit_partition_node[4,1] + circuit_node[4] >= 0
 R1192: - circuit_partition_node[4,2] + circuit_node[4] >= 0
 R1193: - circuit_partition_node[4,3] + circuit_node[4] >= 0
 R1194: - circuit_partition_node[4,4] + circuit_node[4] >= 0
 R1195: - circuit_partition_node[4,5] + circuit_node[4] >= 0
 R1196: - circuit_partition_node[4,6] + circuit_node[4] >= 0
 R1197: - circuit_partition_node[4,7] + circuit_node[4] >= 0
 R1198: - circuit_partition_node[4,8] + circuit_node[4] >= 0
 R1199: - circuit_partition_node[4,9] + circuit_node[4] >= 0
 R1200: - circuit_partition_node[3,0] + circuit_node[3] >= 0
 R1201: - circuit_partition_node[3,1] + circuit_node[3] >= 0
 R1202: - circuit_partition_node[3,2] + circuit_node[3] >= 0
 R1203: - circuit_partition_node[3,3] + circuit_node[3] >= 0
 R1204: - circuit_partition_node[3,4] + circuit_node[3] >= 0
 R1205: - circuit_partition_node[3,5] + circuit_node[3] >= 0
 R1206: - circuit_partition_node[3,6] + circuit_node[3] >= 0
 R1207: - circuit_partition_node[3,7] + circuit_node[3] >= 0
 R1208: - circuit_partition_node[3,8] + circuit_node[3] >= 0
 R1209: - circuit_partition_node[3,9] + circuit_node[3] >= 0
 R1210: - circuit_partition_node[2,0] + circuit_node[2] >= 0
 R1211: - circuit_partition_node[2,1] + circuit_node[2] >= 0
 R1212: - circuit_partition_node[2,2] + circuit_node[2] >= 0
 R1213: - circuit_partition_node[2,3] + circuit_node[2] >= 0
 R1214: - circuit_partition_node[2,4] + circuit_node[2] >= 0
 R1215: - circuit_partition_node[2,5] + circuit_node[2] >= 0
 R1216: - circuit_partition_node[2,6] + circuit_node[2] >= 0
 R1217: - circuit_partition_node[2,7] + circuit_node[2] >= 0
 R1218: - circuit_partition_node[2,8] + circuit_node[2] >= 0
 R1219: - circuit_partition_node[2,9] + circuit_node[2] >= 0
 R1220: - circuit_partition_node[1,0] + circuit_node[1] >= 0
 R1221: - circuit_partition_node[1,1] + circuit_node[1] >= 0
 R1222: - circuit_partition_node[1,2] + circuit_node[1] >= 0
 R1223: - circuit_partition_node[1,3] + circuit_node[1] >= 0
 R1224: - circuit_partition_node[1,4] + circuit_node[1] >= 0
 R1225: - circuit_partition_node[1,5] + circuit_node[1] >= 0
 R1226: - circuit_partition_node[1,6] + circuit_node[1] >= 0
 R1227: - circuit_partition_node[1,7] + circuit_node[1] >= 0
 R1228: - circuit_partition_node[1,8] + circuit_node[1] >= 0
 R1229: - circuit_partition_node[1,9] + circuit_node[1] >= 0
 R1230: - circuit_partition_node[0,0] + circuit_node[0] >= 0
 R1231: - circuit_partition_node[0,1] + circuit_node[0] >= 0
 R1232: - circuit_partition_node[0,2] + circuit_node[0] >= 0
 R1233: - circuit_partition_node[0,3] + circuit_node[0] >= 0
 R1234: - circuit_partition_node[0,4] + circuit_node[0] >= 0
 R1235: - circuit_partition_node[0,5] + circuit_node[0] >= 0
 R1236: - circuit_partition_node[0,6] + circuit_node[0] >= 0
 R1237: - circuit_partition_node[0,7] + circuit_node[0] >= 0
 R1238: - circuit_partition_node[0,8] + circuit_node[0] >= 0
 R1239: - circuit_partition_node[0,9] + circuit_node[0] >= 0
 Circuit_arc_89_-_Both_ends_in_circuit: - circuit_node[5] - circuit_node[4]
   + 2 circuit_arc[4,5] <= 0
 Circuit_arc_88_-_Both_ends_in_circuit: - circuit_node[5] - circuit_node[4]
   + 2 circuit_arc[5,4] <= 0
 Circuit_arc_87_-_Both_ends_in_circuit: - circuit_node[6] - circuit_node[4]
   + 2 circuit_arc[4,6] <= 0
 Circuit_arc_86_-_Both_ends_in_circuit: - circuit_node[6] - circuit_node[4]
   + 2 circuit_arc[6,4] <= 0
 Circuit_arc_85_-_Both_ends_in_circuit: - circuit_node[7] - circuit_node[4]
   + 2 circuit_arc[4,7] <= 0
 Circuit_arc_84_-_Both_ends_in_circuit: - circuit_node[7] - circuit_node[4]
   + 2 circuit_arc[7,4] <= 0
 Circuit_arc_83_-_Both_ends_in_circuit: - circuit_node[8] - circuit_node[4]
   + 2 circuit_arc[4,8] <= 0
 Circuit_arc_82_-_Both_ends_in_circuit: - circuit_node[8] - circuit_node[4]
   + 2 circuit_arc[8,4] <= 0
 Circuit_arc_81_-_Both_ends_in_circuit: - circuit_node[9] - circuit_node[4]
   + 2 circuit_arc[4,9] <= 0
 Circuit_arc_80_-_Both_ends_in_circuit: - circuit_node[9] - circuit_node[4]
   + 2 circuit_arc[9,4] <= 0
 Circuit_arc_79_-_Both_ends_in_circuit: - circuit_node[9] - circuit_node[3]
   + 2 circuit_arc[3,9] <= 0
 Circuit_arc_78_-_Both_ends_in_circuit: - circuit_node[9] - circuit_node[3]
   + 2 circuit_arc[9,3] <= 0
 Circuit_arc_77_-_Both_ends_in_circuit: - circuit_node[8] - circuit_node[3]
   + 2 circuit_arc[3,8] <= 0
 Circuit_arc_76_-_Both_ends_in_circuit: - circuit_node[8] - circuit_node[3]
   + 2 circuit_arc[8,3] <= 0
 Circuit_arc_75_-_Both_ends_in_circuit: - circuit_node[7] - circuit_node[3]
   + 2 circuit_arc[3,7] <= 0
 Circuit_arc_74_-_Both_ends_in_circuit: - circuit_node[7] - circuit_node[3]
   + 2 circuit_arc[7,3] <= 0
 Circuit_arc_73_-_Both_ends_in_circuit: - circuit_node[6] - circuit_node[3]
   + 2 circuit_arc[3,6] <= 0
 Circuit_arc_72_-_Both_ends_in_circuit: - circuit_node[6] - circuit_node[3]
   + 2 circuit_arc[6,3] <= 0
 Circuit_arc_71_-_Both_ends_in_circuit: - circuit_node[5] - circuit_node[3]
   + 2 circuit_arc[3,5] <= 0
 Circuit_arc_70_-_Both_ends_in_circuit: - circuit_node[5] - circuit_node[3]
   + 2 circuit_arc[5,3] <= 0
 Circuit_arc_69_-_Both_ends_in_circuit: - circuit_node[4] - circuit_node[3]
   + 2 circuit_arc[3,4] <= 0
 Circuit_arc_68_-_Both_ends_in_circuit: - circuit_node[4] - circuit_node[3]
   + 2 circuit_arc[4,3] <= 0
 Circuit_arc_67_-_Both_ends_in_circuit: - circuit_node[6] - circuit_node[5]
   + 2 circuit_arc[5,6] <= 0
 Circuit_arc_66_-_Both_ends_in_circuit: - circuit_node[6] - circuit_node[5]
   + 2 circuit_arc[6,5] <= 0
 Circuit_arc_65_-_Both_ends_in_circuit: - circuit_node[7] - circuit_node[5]
   + 2 circuit_arc[5,7] <= 0
 Circuit_arc_64_-_Both_ends_in_circuit: - circuit_node[7] - circuit_node[5]
   + 2 circuit_arc[7,5] <= 0
 Circuit_arc_63_-_Both_ends_in_circuit: - circuit_node[8] - circuit_node[5]
   + 2 circuit_arc[5,8] <= 0
 Circuit_arc_62_-_Both_ends_in_circuit: - circuit_node[8] - circuit_node[5]
   + 2 circuit_arc[8,5] <= 0
 Circuit_arc_61_-_Both_ends_in_circuit: - circuit_node[9] - circuit_node[5]
   + 2 circuit_arc[5,9] <= 0
 Circuit_arc_60_-_Both_ends_in_circuit: - circuit_node[9] - circuit_node[5]
   + 2 circuit_arc[9,5] <= 0
 Circuit_arc_59_-_Both_ends_in_circuit: - circuit_node[9] - circuit_node[2]
   + 2 circuit_arc[2,9] <= 0
 Circuit_arc_58_-_Both_ends_in_circuit: - circuit_node[9] - circuit_node[2]
   + 2 circuit_arc[9,2] <= 0
 Circuit_arc_57_-_Both_ends_in_circuit: - circuit_node[8] - circuit_node[2]
   + 2 circuit_arc[2,8] <= 0
 Circuit_arc_56_-_Both_ends_in_circuit: - circuit_node[8] - circuit_node[2]
   + 2 circuit_arc[8,2] <= 0
 Circuit_arc_55_-_Both_ends_in_circuit: - circuit_node[7] - circuit_node[2]
   + 2 circuit_arc[2,7] <= 0
 Circuit_arc_54_-_Both_ends_in_circuit: - circuit_node[7] - circuit_node[2]
   + 2 circuit_arc[7,2] <= 0
 Circuit_arc_53_-_Both_ends_in_circuit: - circuit_node[6] - circuit_node[2]
   + 2 circuit_arc[2,6] <= 0
 Circuit_arc_52_-_Both_ends_in_circuit: - circuit_node[6] - circuit_node[2]
   + 2 circuit_arc[6,2] <= 0
 Circuit_arc_51_-_Both_ends_in_circuit: - circuit_node[5] - circuit_node[2]
   + 2 circuit_arc[2,5] <= 0
 Circuit_arc_50_-_Both_ends_in_circuit: - circuit_node[5] - circuit_node[2]
   + 2 circuit_arc[5,2] <= 0
 Circuit_arc_49_-_Both_ends_in_circuit: - circuit_node[4] - circuit_node[2]
   + 2 circuit_arc[2,4] <= 0
 Circuit_arc_48_-_Both_ends_in_circuit: - circuit_node[4] - circuit_node[2]
   + 2 circuit_arc[4,2] <= 0
 Circuit_arc_47_-_Both_ends_in_circuit: - circuit_node[3] - circuit_node[2]
   + 2 circuit_arc[2,3] <= 0
 Circuit_arc_46_-_Both_ends_in_circuit: - circuit_node[3] - circuit_node[2]
   + 2 circuit_arc[3,2] <= 0
 Circuit_arc_45_-_Both_ends_in_circuit: - circuit_node[7] - circuit_node[6]
   + 2 circuit_arc[6,7] <= 0
 Circuit_arc_44_-_Both_ends_in_circuit: - circuit_node[7] - circuit_node[6]
   + 2 circuit_arc[7,6] <= 0
 Circuit_arc_43_-_Both_ends_in_circuit: - circuit_node[8] - circuit_node[6]
   + 2 circuit_arc[6,8] <= 0
 Circuit_arc_42_-_Both_ends_in_circuit: - circuit_node[8] - circuit_node[6]
   + 2 circuit_arc[8,6] <= 0
 Circuit_arc_41_-_Both_ends_in_circuit: - circuit_node[9] - circuit_node[6]
   + 2 circuit_arc[6,9] <= 0
 Circuit_arc_40_-_Both_ends_in_circuit: - circuit_node[9] - circuit_node[6]
   + 2 circuit_arc[9,6] <= 0
 Circuit_arc_39_-_Both_ends_in_circuit: - circuit_node[9] - circuit_node[1]
   + 2 circuit_arc[1,9] <= 0
 Circuit_arc_38_-_Both_ends_in_circuit: - circuit_node[9] - circuit_node[1]
   + 2 circuit_arc[9,1] <= 0
 Circuit_arc_37_-_Both_ends_in_circuit: - circuit_node[8] - circuit_node[1]
   + 2 circuit_arc[1,8] <= 0
 Circuit_arc_36_-_Both_ends_in_circuit: - circuit_node[8] - circuit_node[1]
   + 2 circuit_arc[8,1] <= 0
 Circuit_arc_35_-_Both_ends_in_circuit: - circuit_node[7] - circuit_node[1]
   + 2 circuit_arc[1,7] <= 0
 Circuit_arc_34_-_Both_ends_in_circuit: - circuit_node[7] - circuit_node[1]
   + 2 circuit_arc[7,1] <= 0
 Circuit_arc_33_-_Both_ends_in_circuit: - circuit_node[6] - circuit_node[1]
   + 2 circuit_arc[1,6] <= 0
 Circuit_arc_32_-_Both_ends_in_circuit: - circuit_node[6] - circuit_node[1]
   + 2 circuit_arc[6,1] <= 0
 Circuit_arc_31_-_Both_ends_in_circuit: - circuit_node[5] - circuit_node[1]
   + 2 circuit_arc[1,5] <= 0
 Circuit_arc_30_-_Both_ends_in_circuit: - circuit_node[5] - circuit_node[1]
   + 2 circuit_arc[5,1] <= 0
 Circuit_arc_29_-_Both_ends_in_circuit: - circuit_node[4] - circuit_node[1]
   + 2 circuit_arc[1,4] <= 0
 Circuit_arc_28_-_Both_ends_in_circuit: - circuit_node[4] - circuit_node[1]
   + 2 circuit_arc[4,1] <= 0
 Circuit_arc_27_-_Both_ends_in_circuit: - circuit_node[3] - circuit_node[1]
   + 2 circuit_arc[1,3] <= 0
 Circuit_arc_26_-_Both_ends_in_circuit: - circuit_node[3] - circuit_node[1]
   + 2 circuit_arc[3,1] <= 0
 Circuit_arc_25_-_Both_ends_in_circuit: - circuit_node[2] - circuit_node[1]
   + 2 circuit_arc[1,2] <= 0
 Circuit_arc_24_-_Both_ends_in_circuit: - circuit_node[2] - circuit_node[1]
   + 2 circuit_arc[2,1] <= 0
 Circuit_arc_23_-_Both_ends_in_circuit: - circuit_node[8] - circuit_node[7]
   + 2 circuit_arc[7,8] <= 0
 Circuit_arc_22_-_Both_ends_in_circuit: - circuit_node[8] - circuit_node[7]
   + 2 circuit_arc[8,7] <= 0
 Circuit_arc_21_-_Both_ends_in_circuit: - circuit_node[9] - circuit_node[7]
   + 2 circuit_arc[7,9] <= 0
 Circuit_arc_20_-_Both_ends_in_circuit: - circuit_node[9] - circuit_node[7]
   + 2 circuit_arc[9,7] <= 0
 Circuit_arc_19_-_Both_ends_in_circuit: - circuit_node[9] - circuit_node[0]
   + 2 circuit_arc[0,9] <= 0
 Circuit_arc_18_-_Both_ends_in_circuit: - circuit_node[9] - circuit_node[0]
   + 2 circuit_arc[9,0] <= 0
 Circuit_arc_17_-_Both_ends_in_circuit: - circuit_node[8] - circuit_node[0]
   + 2 circuit_arc[0,8] <= 0
 Circuit_arc_16_-_Both_ends_in_circuit: - circuit_node[8] - circuit_node[0]
   + 2 circuit_arc[8,0] <= 0
 Circuit_arc_15_-_Both_ends_in_circuit: - circuit_node[7] - circuit_node[0]
   + 2 circuit_arc[0,7] <= 0
 Circuit_arc_14_-_Both_ends_in_circuit: - circuit_node[7] - circuit_node[0]
   + 2 circuit_arc[7,0] <= 0
 Circuit_arc_13_-_Both_ends_in_circuit: - circuit_node[6] - circuit_node[0]
   + 2 circuit_arc[0,6] <= 0
 Circuit_arc_12_-_Both_ends_in_circuit: - circuit_node[6] - circuit_node[0]
   + 2 circuit_arc[6,0] <= 0
 Circuit_arc_11_-_Both_ends_in_circuit: - circuit_node[5] - circuit_node[0]
   + 2 circuit_arc[0,5] <= 0
 Circuit_arc_10_-_Both_ends_in_circuit: - circuit_node[5] - circuit_node[0]
   + 2 circuit_arc[5,0] <= 0
 Circuit_arc_9_-_Both_ends_in_circuit: - circuit_node[4] - circuit_node[0]
   + 2 circuit_arc[0,4] <= 0
 Circuit_arc_8_-_Both_ends_in_circuit: - circuit_node[4] - circuit_node[0]
   + 2 circuit_arc[4,0] <= 0
 Circuit_arc_7_-_Both_ends_in_circuit: - circuit_node[3] - circuit_node[0]
   + 2 circuit_arc[0,3] <= 0
 Circuit_arc_6_-_Both_ends_in_circuit: - circuit_node[3] - circuit_node[0]
   + 2 circuit_arc[3,0] <= 0
 Circuit_arc_5_-_Both_ends_in_circuit: - circuit_node[2] - circuit_node[0]
   + 2 circuit_arc[0,2] <= 0
 Circuit_arc_4_-_Both_ends_in_circuit: - circuit_node[2] - circuit_node[0]
   + 2 circuit_arc[2,0] <= 0
 Circuit_arc_3_-_Both_ends_in_circuit: - circuit_node[1] - circuit_node[0]
   + 2 circuit_arc[0,1] <= 0
 Circuit_arc_2_-_Both_ends_in_circuit: - circuit_node[1] - circuit_node[0]
   + 2 circuit_arc[1,0] <= 0
 Circuit_arc_1_-_Both_ends_in_circuit: - circuit_node[9] - circuit_node[8]
   + 2 circuit_arc[8,9] <= 0
 Circuit_arc_0_-_Both_ends_in_circuit: - circuit_node[9] - circuit_node[8]
   + 2 circuit_arc[9,8] <= 0
 R1330: circuit_arc[8,9] + circuit_arc[9,8] <= 1
 R1331: circuit_arc[7,9] + circuit_arc[9,7] <= 1
 R1332: circuit_arc[6,9] + circuit_arc[9,6] <= 1
 R1333: circuit_arc[5,9] + circuit_arc[9,5] <= 1
 R1334: circuit_arc[4,9] + circuit_arc[9,4] <= 1
 R1335: circuit_arc[3,9] + circuit_arc[9,3] <= 1
 R1336: circuit_arc[2,9] + circuit_arc[9,2] <= 1
 R1337: circuit_arc[1,9] + circuit_arc[9,1] <= 1
 R1338: circuit_arc[0,9] + circuit_arc[9,0] <= 1
 R1339: circuit_arc[8,9] + circuit_arc[9,8] <= 1
 R1340: circuit_arc[7,8] + circuit_arc[8,7] <= 1
 R1341: circuit_arc[6,8] + circuit_arc[8,6] <= 1
 R1342: circuit_arc[5,8] + circuit_arc[8,5] <= 1
 R1343: circuit_arc[4,8] + circuit_arc[8,4] <= 1
 R1344: circuit_arc[3,8] + circuit_arc[8,3] <= 1
 R1345: circuit_arc[2,8] + circuit_arc[8,2] <= 1
 R1346: circuit_arc[1,8] + circuit_arc[8,1] <= 1
 R1347: circuit_arc[0,8] + circuit_arc[8,0] <= 1
 R1348: circuit_arc[7,9] + circuit_arc[9,7] <= 1
 R1349: circuit_arc[7,8] + circuit_arc[8,7] <= 1
 R1350: circuit_arc[6,7] + circuit_arc[7,6] <= 1
 R1351: circuit_arc[5,7] + circuit_arc[7,5] <= 1
 R1352: circuit_arc[4,7] + circuit_arc[7,4] <= 1
 R1353: circuit_arc[3,7] + circuit_arc[7,3] <= 1
 R1354: circuit_arc[2,7] + circuit_arc[7,2] <= 1
 R1355: circuit_arc[1,7] + circuit_arc[7,1] <= 1
 R1356: circuit_arc[0,7] + circuit_arc[7,0] <= 1
 R1357: circuit_arc[6,9] + circuit_arc[9,6] <= 1
 R1358: circuit_arc[6,8] + circuit_arc[8,6] <= 1
 R1359: circuit_arc[6,7] + circuit_arc[7,6] <= 1
 R1360: circuit_arc[5,6] + circuit_arc[6,5] <= 1
 R1361: circuit_arc[4,6] + circuit_arc[6,4] <= 1
 R1362: circuit_arc[3,6] + circuit_arc[6,3] <= 1
 R1363: circuit_arc[2,6] + circuit_arc[6,2] <= 1
 R1364: circuit_arc[1,6] + circuit_arc[6,1] <= 1
 R1365: circuit_arc[0,6] + circuit_arc[6,0] <= 1
 R1366: circuit_arc[5,9] + circuit_arc[9,5] <= 1
 R1367: circuit_arc[5,8] + circuit_arc[8,5] <= 1
 R1368: circuit_arc[5,7] + circuit_arc[7,5] <= 1
 R1369: circuit_arc[5,6] + circuit_arc[6,5] <= 1
 R1370: circuit_arc[4,5] + circuit_arc[5,4] <= 1
 R1371: circuit_arc[3,5] + circuit_arc[5,3] <= 1
 R1372: circuit_arc[2,5] + circuit_arc[5,2] <= 1
 R1373: circuit_arc[1,5] + circuit_arc[5,1] <= 1
 R1374: circuit_arc[0,5] + circuit_arc[5,0] <= 1
 R1375: circuit_arc[4,9] + circuit_arc[9,4] <= 1
 R1376: circuit_arc[4,8] + circuit_arc[8,4] <= 1
 R1377: circuit_arc[4,7] + circuit_arc[7,4] <= 1
 R1378: circuit_arc[4,6] + circuit_arc[6,4] <= 1
 R1379: circuit_arc[4,5] + circuit_arc[5,4] <= 1
 R1380: circuit_arc[3,4] + circuit_arc[4,3] <= 1
 R1381: circuit_arc[2,4] + circuit_arc[4,2] <= 1
 R1382: circuit_arc[1,4] + circuit_arc[4,1] <= 1
 R1383: circuit_arc[0,4] + circuit_arc[4,0] <= 1
 R1384: circuit_arc[3,9] + circuit_arc[9,3] <= 1
 R1385: circuit_arc[3,8] + circuit_arc[8,3] <= 1
 R1386: circuit_arc[3,7] + circuit_arc[7,3] <= 1
 R1387: circuit_arc[3,6] + circuit_arc[6,3] <= 1
 R1388: circuit_arc[3,5] + circuit_arc[5,3] <= 1
 R1389: circuit_arc[3,4] + circuit_arc[4,3] <= 1
 R1390: circuit_arc[2,3] + circuit_arc[3,2] <= 1
 R1391: circuit_arc[1,3] + circuit_arc[3,1] <= 1
 R1392: circuit_arc[0,3] + circuit_arc[3,0] <= 1
 R1393: circuit_arc[2,9] + circuit_arc[9,2] <= 1
 R1394: circuit_arc[2,8] + circuit_arc[8,2] <= 1
 R1395: circuit_arc[2,7] + circuit_arc[7,2] <= 1
 R1396: circuit_arc[2,6] + circuit_arc[6,2] <= 1
 R1397: circuit_arc[2,5] + circuit_arc[5,2] <= 1
 R1398: circuit_arc[2,4] + circuit_arc[4,2] <= 1
 R1399: circuit_arc[2,3] + circuit_arc[3,2] <= 1
 R1400: circuit_arc[1,2] + circuit_arc[2,1] <= 1
 R1401: circuit_arc[0,2] + circuit_arc[2,0] <= 1
 R1402: circuit_arc[1,9] + circuit_arc[9,1] <= 1
 R1403: circuit_arc[1,8] + circuit_arc[8,1] <= 1
 R1404: circuit_arc[1,7] + circuit_arc[7,1] <= 1
 R1405: circuit_arc[1,6] + circuit_arc[6,1] <= 1
 R1406: circuit_arc[1,5] + circuit_arc[5,1] <= 1
 R1407: circuit_arc[1,4] + circuit_arc[4,1] <= 1
 R1408: circuit_arc[1,3] + circuit_arc[3,1] <= 1
 R1409: circuit_arc[1,2] + circuit_arc[2,1] <= 1
 R1410: circuit_arc[0,1] + circuit_arc[1,0] <= 1
 R1411: circuit_arc[0,9] + circuit_arc[9,0] <= 1
 R1412: circuit_arc[0,8] + circuit_arc[8,0] <= 1
 R1413: circuit_arc[0,7] + circuit_arc[7,0] <= 1
 R1414: circuit_arc[0,6] + circuit_arc[6,0] <= 1
 R1415: circuit_arc[0,5] + circuit_arc[5,0] <= 1
 R1416: circuit_arc[0,4] + circuit_arc[4,0] <= 1
 R1417: circuit_arc[0,3] + circuit_arc[3,0] <= 1
 R1418: circuit_arc[0,2] + circuit_arc[2,0] <= 1
 R1419: circuit_arc[0,1] + circuit_arc[1,0] <= 1
 Circuit_node_9_-_In-Degree_1: - circuit_node[9] + circuit_arc[4,9]
   + circuit_arc[3,9] + circuit_arc[5,9] + circuit_arc[2,9]
   + circuit_arc[6,9] + circuit_arc[1,9] + circuit_arc[7,9]
   + circuit_arc[0,9] + circuit_arc[8,9] = 0
 Circuit_node_9_-_Out-Degree_1: - circuit_node[9] + circuit_arc[9,4]
   + circuit_arc[9,3] + circuit_arc[9,5] + circuit_arc[9,2]
   + circuit_arc[9,6] + circuit_arc[9,1] + circuit_arc[9,7]
   + circuit_arc[9,0] + circuit_arc[9,8] = 0
 Circuit_node_8_-_In-Degree_1: - circuit_node[8] + circuit_arc[4,8]
   + circuit_arc[3,8] + circuit_arc[5,8] + circuit_arc[2,8]
   + circuit_arc[6,8] + circuit_arc[1,8] + circuit_arc[7,8]
   + circuit_arc[0,8] + circuit_arc[9,8] = 0
 Circuit_node_8_-_Out-Degree_1: - circuit_node[8] + circuit_arc[8,4]
   + circuit_arc[8,3] + circuit_arc[8,5] + circuit_arc[8,2]
   + circuit_arc[8,6] + circuit_arc[8,1] + circuit_arc[8,7]
   + circuit_arc[8,0] + circuit_arc[8,9] = 0
 Circuit_node_7_-_In-Degree_1: - circuit_node[7] + circuit_arc[4,7]
   + circuit_arc[3,7] + circuit_arc[5,7] + circuit_arc[2,7]
   + circuit_arc[6,7] + circuit_arc[1,7] + circuit_arc[8,7]
   + circuit_arc[9,7] + circuit_arc[0,7] = 0
 Circuit_node_7_-_Out-Degree_1: - circuit_node[7] + circuit_arc[7,4]
   + circuit_arc[7,3] + circuit_arc[7,5] + circuit_arc[7,2]
   + circuit_arc[7,6] + circuit_arc[7,1] + circuit_arc[7,8]
   + circuit_arc[7,9] + circuit_arc[7,0] = 0
 Circuit_node_6_-_In-Degree_1: - circuit_node[6] + circuit_arc[4,6]
   + circuit_arc[3,6] + circuit_arc[5,6] + circuit_arc[2,6]
   + circuit_arc[7,6] + circuit_arc[8,6] + circuit_arc[9,6]
   + circuit_arc[1,6] + circuit_arc[0,6] = 0
 Circuit_node_6_-_Out-Degree_1: - circuit_node[6] + circuit_arc[6,4]
   + circuit_arc[6,3] + circuit_arc[6,5] + circuit_arc[6,2]
   + circuit_arc[6,7] + circuit_arc[6,8] + circuit_arc[6,9]
   + circuit_arc[6,1] + circuit_arc[6,0] = 0
 Circuit_node_5_-_In-Degree_1: - circuit_node[5] + circuit_arc[4,5]
   + circuit_arc[3,5] + circuit_arc[6,5] + circuit_arc[7,5]
   + circuit_arc[8,5] + circuit_arc[9,5] + circuit_arc[2,5]
   + circuit_arc[1,5] + circuit_arc[0,5] = 0
 Circuit_node_5_-_Out-Degree_1: - circuit_node[5] + circuit_arc[5,4]
   + circuit_arc[5,3] + circuit_arc[5,6] + circuit_arc[5,7]
   + circuit_arc[5,8] + circuit_arc[5,9] + circuit_arc[5,2]
   + circuit_arc[5,1] + circuit_arc[5,0] = 0
 Circuit_node_4_-_In-Degree_1: - circuit_node[4] + circuit_arc[5,4]
   + circuit_arc[6,4] + circuit_arc[7,4] + circuit_arc[8,4]
   + circuit_arc[9,4] + circuit_arc[3,4] + circuit_arc[2,4]
   + circuit_arc[1,4] + circuit_arc[0,4] = 0
 Circuit_node_4_-_Out-Degree_1: - circuit_node[4] + circuit_arc[4,5]
   + circuit_arc[4,6] + circuit_arc[4,7] + circuit_arc[4,8]
   + circuit_arc[4,9] + circuit_arc[4,3] + circuit_arc[4,2]
   + circuit_arc[4,1] + circuit_arc[4,0] = 0
 Circuit_node_3_-_In-Degree_1: - circuit_node[3] + circuit_arc[9,3]
   + circuit_arc[8,3] + circuit_arc[7,3] + circuit_arc[6,3]
   + circuit_arc[5,3] + circuit_arc[4,3] + circuit_arc[2,3]
   + circuit_arc[1,3] + circuit_arc[0,3] = 0
 Circuit_node_3_-_Out-Degree_1: - circuit_node[3] + circuit_arc[3,9]
   + circuit_arc[3,8] + circuit_arc[3,7] + circuit_arc[3,6]
   + circuit_arc[3,5] + circuit_arc[3,4] + circuit_arc[3,2]
   + circuit_arc[3,1] + circuit_arc[3,0] = 0
 Circuit_node_2_-_In-Degree_1: - circuit_node[2] + circuit_arc[9,2]
   + circuit_arc[8,2] + circuit_arc[7,2] + circuit_arc[6,2]
   + circuit_arc[5,2] + circuit_arc[4,2] + circuit_arc[3,2]
   + circuit_arc[1,2] + circuit_arc[0,2] = 0
 Circuit_node_2_-_Out-Degree_1: - circuit_node[2] + circuit_arc[2,9]
   + circuit_arc[2,8] + circuit_arc[2,7] + circuit_arc[2,6]
   + circuit_arc[2,5] + circuit_arc[2,4] + circuit_arc[2,3]
   + circuit_arc[2,1] + circuit_arc[2,0] = 0
 Circuit_node_1_-_In-Degree_1: - circuit_node[1] + circuit_arc[9,1]
   + circuit_arc[8,1] + circuit_arc[7,1] + circuit_arc[6,1]
   + circuit_arc[5,1] + circuit_arc[4,1] + circuit_arc[3,1]
   + circuit_arc[2,1] + circuit_arc[0,1] = 0
 Circuit_node_1_-_Out-Degree_1: - circuit_node[1] + circuit_arc[1,9]
   + circuit_arc[1,8] + circuit_arc[1,7] + circuit_arc[1,6]
   + circuit_arc[1,5] + circuit_arc[1,4] + circuit_arc[1,3]
   + circuit_arc[1,2] + circuit_arc[1,0] = 0
 Circuit_node_0_-_In-Degree_1: - circuit_node[0] + circuit_arc[9,0]
   + circuit_arc[8,0] + circuit_arc[7,0] + circuit_arc[6,0]
   + circuit_arc[5,0] + circuit_arc[4,0] + circuit_arc[3,0]
   + circuit_arc[2,0] + circuit_arc[1,0] = 0
 Circuit_node_0_-_Out-Degree_1: - circuit_node[0] + circuit_arc[0,9]
   + circuit_arc[0,8] + circuit_arc[0,7] + circuit_arc[0,6]
   + circuit_arc[0,5] + circuit_arc[0,4] + circuit_arc[0,3]
   + circuit_arc[0,2] + circuit_arc[0,1] = 0
 R1440: - partition_used[0] - partition_used[1] - partition_used[2]
   - partition_used[3] - partition_used[4] - partition_used[5]
   - partition_used[6] - partition_used[7] - partition_used[8]
   - partition_used[9] + circuit_order[9] <= -1
 R1441: - partition_used[0] - partition_used[1] - partition_used[2]
   - partition_used[3] - partition_used[4] - partition_used[5]
   - partition_used[6] - partition_used[7] - partition_used[8]
   - partition_used[9] + circuit_order[8] <= -1
 R1442: - partition_used[0] - partition_used[1] - partition_used[2]
   - partition_used[3] - partition_used[4] - partition_used[5]
   - partition_used[6] - partition_used[7] - partition_used[8]
   - partition_used[9] + circuit_order[7] <= -1
 R1443: - partition_used[0] - partition_used[1] - partition_used[2]
   - partition_used[3] - partition_used[4] - partition_used[5]
   - partition_used[6] - partition_used[7] - partition_used[8]
   - partition_used[9] + circuit_order[6] <= -1
 R1444: - partition_used[0] - partition_used[1] - partition_used[2]
   - partition_used[3] - partition_used[4] - partition_used[5]
   - partition_used[6] - partition_used[7] - partition_used[8]
   - partition_used[9] + circuit_order[5] <= -1
 R1445: - partition_used[0] - partition_used[1] - partition_used[2]
   - partition_used[3] - partition_used[4] - partition_used[5]
   - partition_used[6] - partition_used[7] - partition_used[8]
   - partition_used[9] + circuit_order[4] <= -1
 R1446: - partition_used[0] - partition_used[1] - partition_used[2]
   - partition_used[3] - partition_used[4] - partition_used[5]
   - partition_used[6] - partition_used[7] - partition_used[8]
   - partition_used[9] + circuit_order[3] <= -1
 R1447: - partition_used[0] - partition_used[1] - partition_used[2]
   - partition_used[3] - partition_used[4] - partition_used[5]
   - partition_used[6] - partition_used[7] - partition_used[8]
   - partition_used[9] + circuit_order[2] <= -1
 R1448: - partition_used[0] - partition_used[1] - partition_used[2]
   - partition_used[3] - partition_used[4] - partition_used[5]
   - partition_used[6] - partition_used[7] - partition_used[8]
   - partition_used[9] + circuit_order[1] <= -1
 R1449: - partition_used[0] - partition_used[1] - partition_used[2]
   - partition_used[3] - partition_used[4] - partition_used[5]
   - partition_used[6] - partition_used[7] - partition_used[8]
   - partition_used[9] + circuit_order[0] <= -1
 Circuit_arc_89_-_Order: 10 circuit_partition_node[5,0] + circuit_order[5]
   - circuit_order[4] - 10 circuit_arc[4,5] >= -9
 Circuit_arc_88_-_Order: 10 circuit_partition_node[4,0] - circuit_order[5]
   + circuit_order[4] - 10 circuit_arc[5,4] >= -9
 Circuit_arc_87_-_Order: 10 circuit_partition_node[6,0] + circuit_order[6]
   - circuit_order[4] - 10 circuit_arc[4,6] >= -9
 Circuit_arc_86_-_Order: 10 circuit_partition_node[4,0] - circuit_order[6]
   + circuit_order[4] - 10 circuit_arc[6,4] >= -9
 Circuit_arc_85_-_Order: 10 circuit_partition_node[7,0] + circuit_order[7]
   - circuit_order[4] - 10 circuit_arc[4,7] >= -9
 Circuit_arc_84_-_Order: 10 circuit_partition_node[4,0] - circuit_order[7]
   + circuit_order[4] - 10 circuit_arc[7,4] >= -9
 Circuit_arc_83_-_Order: 10 circuit_partition_node[8,0] + circuit_order[8]
   - circuit_order[4] - 10 circuit_arc[4,8] >= -9
 Circuit_arc_82_-_Order: 10 circuit_partition_node[4,0] - circuit_order[8]
   + circuit_order[4] - 10 circuit_arc[8,4] >= -9
 Circuit_arc_81_-_Order: 10 circuit_partition_node[9,0] + circuit_order[9]
   - circuit_order[4] - 10 circuit_arc[4,9] >= -9
 Circuit_arc_80_-_Order: 10 circuit_partition_node[4,0] - circuit_order[9]
   + circuit_order[4] - 10 circuit_arc[9,4] >= -9
 Circuit_arc_79_-_Order: 10 circuit_partition_node[9,0] + circuit_order[9]
   - circuit_order[3] - 10 circuit_arc[3,9] >= -9
 Circuit_arc_78_-_Order: 10 circuit_partition_node[3,0] - circuit_order[9]
   + circuit_order[3] - 10 circuit_arc[9,3] >= -9
 Circuit_arc_77_-_Order: 10 circuit_partition_node[8,0] + circuit_order[8]
   - circuit_order[3] - 10 circuit_arc[3,8] >= -9
 Circuit_arc_76_-_Order: 10 circuit_partition_node[3,0] - circuit_order[8]
   + circuit_order[3] - 10 circuit_arc[8,3] >= -9
 Circuit_arc_75_-_Order: 10 circuit_partition_node[7,0] + circuit_order[7]
   - circuit_order[3] - 10 circuit_arc[3,7] >= -9
 Circuit_arc_74_-_Order: 10 circuit_partition_node[3,0] - circuit_order[7]
   + circuit_order[3] - 10 circuit_arc[7,3] >= -9
 Circuit_arc_73_-_Order: 10 circuit_partition_node[6,0] + circuit_order[6]
   - circuit_order[3] - 10 circuit_arc[3,6] >= -9
 Circuit_arc_72_-_Order: 10 circuit_partition_node[3,0] - circuit_order[6]
   + circuit_order[3] - 10 circuit_arc[6,3] >= -9
 Circuit_arc_71_-_Order: 10 circuit_partition_node[5,0] + circuit_order[5]
   - circuit_order[3] - 10 circuit_arc[3,5] >= -9
 Circuit_arc_70_-_Order: 10 circuit_partition_node[3,0] - circuit_order[5]
   + circuit_order[3] - 10 circuit_arc[5,3] >= -9
 Circuit_arc_69_-_Order: 10 circuit_partition_node[4,0] + circuit_order[4]
   - circuit_order[3] - 10 circuit_arc[3,4] >= -9
 Circuit_arc_68_-_Order: 10 circuit_partition_node[3,0] - circuit_order[4]
   + circuit_order[3] - 10 circuit_arc[4,3] >= -9
 Circuit_arc_67_-_Order: 10 circuit_partition_node[6,0] + circuit_order[6]
   - circuit_order[5] - 10 circuit_arc[5,6] >= -9
 Circuit_arc_66_-_Order: 10 circuit_partition_node[5,0] - circuit_order[6]
   + circuit_order[5] - 10 circuit_arc[6,5] >= -9
 Circuit_arc_65_-_Order: 10 circuit_partition_node[7,0] + circuit_order[7]
   - circuit_order[5] - 10 circuit_arc[5,7] >= -9
 Circuit_arc_64_-_Order: 10 circuit_partition_node[5,0] - circuit_order[7]
   + circuit_order[5] - 10 circuit_arc[7,5] >= -9
 Circuit_arc_63_-_Order: 10 circuit_partition_node[8,0] + circuit_order[8]
   - circuit_order[5] - 10 circuit_arc[5,8] >= -9
 Circuit_arc_62_-_Order: 10 circuit_partition_node[5,0] - circuit_order[8]
   + circuit_order[5] - 10 circuit_arc[8,5] >= -9
 Circuit_arc_61_-_Order: 10 circuit_partition_node[9,0] + circuit_order[9]
   - circuit_order[5] - 10 circuit_arc[5,9] >= -9
 Circuit_arc_60_-_Order: 10 circuit_partition_node[5,0] - circuit_order[9]
   + circuit_order[5] - 10 circuit_arc[9,5] >= -9
 Circuit_arc_59_-_Order: 10 circuit_partition_node[9,0] + circuit_order[9]
   - circuit_order[2] - 10 circuit_arc[2,9] >= -9
 Circuit_arc_58_-_Order: 10 circuit_partition_node[2,0] - circuit_order[9]
   + circuit_order[2] - 10 circuit_arc[9,2] >= -9
 Circuit_arc_57_-_Order: 10 circuit_partition_node[8,0] + circuit_order[8]
   - circuit_order[2] - 10 circuit_arc[2,8] >= -9
 Circuit_arc_56_-_Order: 10 circuit_partition_node[2,0] - circuit_order[8]
   + circuit_order[2] - 10 circuit_arc[8,2] >= -9
 Circuit_arc_55_-_Order: 10 circuit_partition_node[7,0] + circuit_order[7]
   - circuit_order[2] - 10 circuit_arc[2,7] >= -9
 Circuit_arc_54_-_Order: 10 circuit_partition_node[2,0] - circuit_order[7]
   + circuit_order[2] - 10 circuit_arc[7,2] >= -9
 Circuit_arc_53_-_Order: 10 circuit_partition_node[6,0] + circuit_order[6]
   - circuit_order[2] - 10 circuit_arc[2,6] >= -9
 Circuit_arc_52_-_Order: 10 circuit_partition_node[2,0] - circuit_order[6]
   + circuit_order[2] - 10 circuit_arc[6,2] >= -9
 Circuit_arc_51_-_Order: 10 circuit_partition_node[5,0] + circuit_order[5]
   - circuit_order[2] - 10 circuit_arc[2,5] >= -9
 Circuit_arc_50_-_Order: 10 circuit_partition_node[2,0] - circuit_order[5]
   + circuit_order[2] - 10 circuit_arc[5,2] >= -9
 Circuit_arc_49_-_Order: 10 circuit_partition_node[4,0] + circuit_order[4]
   - circuit_order[2] - 10 circuit_arc[2,4] >= -9
 Circuit_arc_48_-_Order: 10 circuit_partition_node[2,0] - circuit_order[4]
   + circuit_order[2] - 10 circuit_arc[4,2] >= -9
 Circuit_arc_47_-_Order: 10 circuit_partition_node[3,0] + circuit_order[3]
   - circuit_order[2] - 10 circuit_arc[2,3] >= -9
 Circuit_arc_46_-_Order: 10 circuit_partition_node[2,0] - circuit_order[3]
   + circuit_order[2] - 10 circuit_arc[3,2] >= -9
 Circuit_arc_45_-_Order: 10 circuit_partition_node[7,0] + circuit_order[7]
   - circuit_order[6] - 10 circuit_arc[6,7] >= -9
 Circuit_arc_44_-_Order: 10 circuit_partition_node[6,0] - circuit_order[7]
   + circuit_order[6] - 10 circuit_arc[7,6] >= -9
 Circuit_arc_43_-_Order: 10 circuit_partition_node[8,0] + circuit_order[8]
   - circuit_order[6] - 10 circuit_arc[6,8] >= -9
 Circuit_arc_42_-_Order: 10 circuit_partition_node[6,0] - circuit_order[8]
   + circuit_order[6] - 10 circuit_arc[8,6] >= -9
 Circuit_arc_41_-_Order: 10 circuit_partition_node[9,0] + circuit_order[9]
   - circuit_order[6] - 10 circuit_arc[6,9] >= -9
 Circuit_arc_40_-_Order: 10 circuit_partition_node[6,0] - circuit_order[9]
   + circuit_order[6] - 10 circuit_arc[9,6] >= -9
 Circuit_arc_39_-_Order: 10 circuit_partition_node[9,0] + circuit_order[9]
   - circuit_order[1] - 10 circuit_arc[1,9] >= -9
 Circuit_arc_38_-_Order: 10 circuit_partition_node[1,0] - circuit_order[9]
   + circuit_order[1] - 10 circuit_arc[9,1] >= -9
 Circuit_arc_37_-_Order: 10 circuit_partition_node[8,0] + circuit_order[8]
   - circuit_order[1] - 10 circuit_arc[1,8] >= -9
 Circuit_arc_36_-_Order: 10 circuit_partition_node[1,0] - circuit_order[8]
   + circuit_order[1] - 10 circuit_arc[8,1] >= -9
 Circuit_arc_35_-_Order: 10 circuit_partition_node[7,0] + circuit_order[7]
   - circuit_order[1] - 10 circuit_arc[1,7] >= -9
 Circuit_arc_34_-_Order: 10 circuit_partition_node[1,0] - circuit_order[7]
   + circuit_order[1] - 10 circuit_arc[7,1] >= -9
 Circuit_arc_33_-_Order: 10 circuit_partition_node[6,0] + circuit_order[6]
   - circuit_order[1] - 10 circuit_arc[1,6] >= -9
 Circuit_arc_32_-_Order: 10 circuit_partition_node[1,0] - circuit_order[6]
   + circuit_order[1] - 10 circuit_arc[6,1] >= -9
 Circuit_arc_31_-_Order: 10 circuit_partition_node[5,0] + circuit_order[5]
   - circuit_order[1] - 10 circuit_arc[1,5] >= -9
 Circuit_arc_30_-_Order: 10 circuit_partition_node[1,0] - circuit_order[5]
   + circuit_order[1] - 10 circuit_arc[5,1] >= -9
 Circuit_arc_29_-_Order: 10 circuit_partition_node[4,0] + circuit_order[4]
   - circuit_order[1] - 10 circuit_arc[1,4] >= -9
 Circuit_arc_28_-_Order: 10 circuit_partition_node[1,0] - circuit_order[4]
   + circuit_order[1] - 10 circuit_arc[4,1] >= -9
 Circuit_arc_27_-_Order: 10 circuit_partition_node[3,0] + circuit_order[3]
   - circuit_order[1] - 10 circuit_arc[1,3] >= -9
 Circuit_arc_26_-_Order: 10 circuit_partition_node[1,0] - circuit_order[3]
   + circuit_order[1] - 10 circuit_arc[3,1] >= -9
 Circuit_arc_25_-_Order: 10 circuit_partition_node[2,0] + circuit_order[2]
   - circuit_order[1] - 10 circuit_arc[1,2] >= -9
 Circuit_arc_24_-_Order: 10 circuit_partition_node[1,0] - circuit_order[2]
   + circuit_order[1] - 10 circuit_arc[2,1] >= -9
 Circuit_arc_23_-_Order: 10 circuit_partition_node[8,0] + circuit_order[8]
   - circuit_order[7] - 10 circuit_arc[7,8] >= -9
 Circuit_arc_22_-_Order: 10 circuit_partition_node[7,0] - circuit_order[8]
   + circuit_order[7] - 10 circuit_arc[8,7] >= -9
 Circuit_arc_21_-_Order: 10 circuit_partition_node[9,0] + circuit_order[9]
   - circuit_order[7] - 10 circuit_arc[7,9] >= -9
 Circuit_arc_20_-_Order: 10 circuit_partition_node[7,0] - circuit_order[9]
   + circuit_order[7] - 10 circuit_arc[9,7] >= -9
 Circuit_arc_19_-_Order: 10 circuit_partition_node[9,0] + circuit_order[9]
   - circuit_order[0] - 10 circuit_arc[0,9] >= -9
 Circuit_arc_18_-_Order: 10 circuit_partition_node[0,0] - circuit_order[9]
   + circuit_order[0] - 10 circuit_arc[9,0] >= -9
 Circuit_arc_17_-_Order: 10 circuit_partition_node[8,0] + circuit_order[8]
   - circuit_order[0] - 10 circuit_arc[0,8] >= -9
 Circuit_arc_16_-_Order: 10 circuit_partition_node[0,0] - circuit_order[8]
   + circuit_order[0] - 10 circuit_arc[8,0] >= -9
 Circuit_arc_15_-_Order: 10 circuit_partition_node[7,0] + circuit_order[7]
   - circuit_order[0] - 10 circuit_arc[0,7] >= -9
 Circuit_arc_14_-_Order: 10 circuit_partition_node[0,0] - circuit_order[7]
   + circuit_order[0] - 10 circuit_arc[7,0] >= -9
 Circuit_arc_13_-_Order: 10 circuit_partition_node[6,0] + circuit_order[6]
   - circuit_order[0] - 10 circuit_arc[0,6] >= -9
 Circuit_arc_12_-_Order: 10 circuit_partition_node[0,0] - circuit_order[6]
   + circuit_order[0] - 10 circuit_arc[6,0] >= -9
 Circuit_arc_11_-_Order: 10 circuit_partition_node[5,0] + circuit_order[5]
   - circuit_order[0] - 10 circuit_arc[0,5] >= -9
 Circuit_arc_10_-_Order: 10 circuit_partition_node[0,0] - circuit_order[5]
   + circuit_order[0] - 10 circuit_arc[5,0] >= -9
 Circuit_arc_9_-_Order: 10 circuit_partition_node[4,0] + circuit_order[4]
   - circuit_order[0] - 10 circuit_arc[0,4] >= -9
 Circuit_arc_8_-_Order: 10 circuit_partition_node[0,0] - circuit_order[4]
   + circuit_order[0] - 10 circuit_arc[4,0] >= -9
 Circuit_arc_7_-_Order: 10 circuit_partition_node[3,0] + circuit_order[3]
   - circuit_order[0] - 10 circuit_arc[0,3] >= -9
 Circuit_arc_6_-_Order: 10 circuit_partition_node[0,0] - circuit_order[3]
   + circuit_order[0] - 10 circuit_arc[3,0] >= -9
 Circuit_arc_5_-_Order: 10 circuit_partition_node[2,0] + circuit_order[2]
   - circuit_order[0] - 10 circuit_arc[0,2] >= -9
 Circuit_arc_4_-_Order: 10 circuit_partition_node[0,0] - circuit_order[2]
   + circuit_order[0] - 10 circuit_arc[2,0] >= -9
 Circuit_arc_3_-_Order: 10 circuit_partition_node[1,0] + circuit_order[1]
   - circuit_order[0] - 10 circuit_arc[0,1] >= -9
 Circuit_arc_2_-_Order: 10 circuit_partition_node[0,0] - circuit_order[1]
   + circuit_order[0] - 10 circuit_arc[1,0] >= -9
 Circuit_arc_1_-_Order: 10 circuit_partition_node[9,0] + circuit_order[9]
   - circuit_order[8] - 10 circuit_arc[8,9] >= -9
 Circuit_arc_0_-_Order: 10 circuit_partition_node[8,0] - circuit_order[9]
   + circuit_order[8] - 10 circuit_arc[9,8] >= -9
 Circuit_arc_0_-_Order: - circuit_partition_node[5,0]
   + circuit_partition_node[4,0] - circuit_node[4] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[6,0]
   + circuit_partition_node[4,0] - circuit_node[4] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[7,0]
   + circuit_partition_node[4,0] - circuit_node[4] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[8,0]
   + circuit_partition_node[4,0] - circuit_node[4] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[9,0]
   + circuit_partition_node[4,0] - circuit_node[4] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[9,0]
   + circuit_partition_node[3,0] - circuit_node[3] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[8,0]
   + circuit_partition_node[3,0] - circuit_node[3] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[7,0]
   + circuit_partition_node[3,0] - circuit_node[3] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[6,0]
   + circuit_partition_node[3,0] - circuit_node[3] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[5,0]
   + circuit_partition_node[3,0] - circuit_node[3] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[4,0]
   + circuit_partition_node[3,0] - circuit_node[3] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[6,0]
   + circuit_partition_node[5,0] - circuit_node[5] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[7,0]
   + circuit_partition_node[5,0] - circuit_node[5] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[8,0]
   + circuit_partition_node[5,0] - circuit_node[5] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[9,0]
   + circuit_partition_node[5,0] - circuit_node[5] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[9,0]
   + circuit_partition_node[2,0] - circuit_node[2] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[8,0]
   + circuit_partition_node[2,0] - circuit_node[2] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[7,0]
   + circuit_partition_node[2,0] - circuit_node[2] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[6,0]
   + circuit_partition_node[2,0] - circuit_node[2] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[5,0]
   + circuit_partition_node[2,0] - circuit_node[2] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[4,0]
   + circuit_partition_node[2,0] - circuit_node[2] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[3,0]
   + circuit_partition_node[2,0] - circuit_node[2] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[7,0]
   + circuit_partition_node[6,0] - circuit_node[6] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[8,0]
   + circuit_partition_node[6,0] - circuit_node[6] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[9,0]
   + circuit_partition_node[6,0] - circuit_node[6] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[9,0]
   + circuit_partition_node[1,0] - circuit_node[1] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[8,0]
   + circuit_partition_node[1,0] - circuit_node[1] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[7,0]
   + circuit_partition_node[1,0] - circuit_node[1] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[6,0]
   + circuit_partition_node[1,0] - circuit_node[1] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[5,0]
   + circuit_partition_node[1,0] - circuit_node[1] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[4,0]
   + circuit_partition_node[1,0] - circuit_node[1] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[3,0]
   + circuit_partition_node[1,0] - circuit_node[1] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[2,0]
   + circuit_partition_node[1,0] - circuit_node[1] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[8,0]
   + circuit_partition_node[7,0] - circuit_node[7] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[9,0]
   + circuit_partition_node[7,0] - circuit_node[7] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[9,0]
   + circuit_partition_node[0,0] - circuit_node[0] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[8,0]
   + circuit_partition_node[0,0] - circuit_node[0] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[7,0]
   + circuit_partition_node[0,0] - circuit_node[0] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[6,0]
   + circuit_partition_node[0,0] - circuit_node[0] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[5,0]
   + circuit_partition_node[0,0] - circuit_node[0] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[4,0]
   + circuit_partition_node[0,0] - circuit_node[0] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[3,0]
   + circuit_partition_node[0,0] - circuit_node[0] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[2,0]
   + circuit_partition_node[0,0] - circuit_node[0] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[1,0]
   + circuit_partition_node[0,0] - circuit_node[0] >= -1
 Circuit_arc_0_-_Order: - circuit_partition_node[9,0]
   + circuit_partition_node[8,0] - circuit_node[8] >= -1
Bounds
 circuit_order[9] <= 9
 circuit_order[8] <= 9
 circuit_order[7] <= 9
 circuit_order[6] <= 9
 circuit_order[5] <= 9
 circuit_order[4] <= 9
 circuit_order[3] <= 9
 circuit_order[2] <= 9
 circuit_order[1] <= 9
 circuit_order[0] <= 9
Binaries
 partition[9,0] partition[9,1] partition[9,2] partition[9,3] partition[9,4]
 partition[9,5] partition[9,6] partition[9,7] partition[9,8] partition[9,9]
 partition[8,0] partition[8,1] partition[8,2] partition[8,3] partition[8,4]
 partition[8,5] partition[8,6] partition[8,7] partition[8,8] partition[8,9]
 partition[7,0] partition[7,1] partition[7,2] partition[7,3] partition[7,4]
 partition[7,5] partition[7,6] partition[7,7] partition[7,8] partition[7,9]
 partition[6,0] partition[6,1] partition[6,2] partition[6,3] partition[6,4]
 partition[6,5] partition[6,6] partition[6,7] partition[6,8] partition[6,9]
 partition[5,0] partition[5,1] partition[5,2] partition[5,3] partition[5,4]
 partition[5,5] partition[5,6] partition[5,7] partition[5,8] partition[5,9]
 partition[4,0] partition[4,1] partition[4,2] partition[4,3] partition[4,4]
 partition[4,5] partition[4,6] partition[4,7] partition[4,8] partition[4,9]
 partition[3,0] partition[3,1] partition[3,2] partition[3,3] partition[3,4]
 partition[3,5] partition[3,6] partition[3,7] partition[3,8] partition[3,9]
 partition[2,0] partition[2,1] partition[2,2] partition[2,3] partition[2,4]
 partition[2,5] partition[2,6] partition[2,7] partition[2,8] partition[2,9]
 partition[1,0] partition[1,1] partition[1,2] partition[1,3] partition[1,4]
 partition[1,5] partition[1,6] partition[1,7] partition[1,8] partition[1,9]
 partition[0,0] partition[0,1] partition[0,2] partition[0,3] partition[0,4]
 partition[0,5] partition[0,6] partition[0,7] partition[0,8] partition[0,9]
 partition_used[0] partition_used[1] partition_used[2] partition_used[3]
 partition_used[4] partition_used[5] partition_used[6] partition_used[7]
 partition_used[8] partition_used[9] circuit_partition_node[9,0]
 circuit_partition_node[9,1] circuit_partition_node[9,2]
 circuit_partition_node[9,3] circuit_partition_node[9,4]
 circuit_partition_node[9,5] circuit_partition_node[9,6]
 circuit_partition_node[9,7] circuit_partition_node[9,8]
 circuit_partition_node[9,9] circuit_node[9] circuit_partition_node[8,0]
 circuit_partition_node[8,1] circuit_partition_node[8,2]
 circuit_partition_node[8,3] circuit_partition_node[8,4]
 circuit_partition_node[8,5] circuit_partition_node[8,6]
 circuit_partition_node[8,7] circuit_partition_node[8,8]
 circuit_partition_node[8,9] circuit_node[8] circuit_partition_node[7,0]
 circuit_partition_node[7,1] circuit_partition_node[7,2]
 circuit_partition_node[7,3] circuit_partition_node[7,4]
 circuit_partition_node[7,5] circuit_partition_node[7,6]
 circuit_partition_node[7,7] circuit_partition_node[7,8]
 circuit_partition_node[7,9] circuit_node[7] circuit_partition_node[6,0]
 circuit_partition_node[6,1] circuit_partition_node[6,2]
 circuit_partition_node[6,3] circuit_partition_node[6,4]
 circuit_partition_node[6,5] circuit_partition_node[6,6]
 circuit_partition_node[6,7] circuit_partition_node[6,8]
 circuit_partition_node[6,9] circuit_node[6] circuit_partition_node[5,0]
 circuit_partition_node[5,1] circuit_partition_node[5,2]
 circuit_partition_node[5,3] circuit_partition_node[5,4]
 circuit_partition_node[5,5] circuit_partition_node[5,6]
 circuit_partition_node[5,7] circuit_partition_node[5,8]
 circuit_partition_node[5,9] circuit_node[5] circuit_partition_node[4,0]
 circuit_partition_node[4,1] circuit_partition_node[4,2]
 circuit_partition_node[4,3] circuit_partition_node[4,4]
 circuit_partition_node[4,5] circuit_partition_node[4,6]
 circuit_partition_node[4,7] circuit_partition_node[4,8]
 circuit_partition_node[4,9] circuit_node[4] circuit_partition_node[3,0]
 circuit_partition_node[3,1] circuit_partition_node[3,2]
 circuit_partition_node[3,3] circuit_partition_node[3,4]
 circuit_partition_node[3,5] circuit_partition_node[3,6]
 circuit_partition_node[3,7] circuit_partition_node[3,8]
 circuit_partition_node[3,9] circuit_node[3] circuit_partition_node[2,0]
 circuit_partition_node[2,1] circuit_partition_node[2,2]
 circuit_partition_node[2,3] circuit_partition_node[2,4]
 circuit_partition_node[2,5] circuit_partition_node[2,6]
 circuit_partition_node[2,7] circuit_partition_node[2,8]
 circuit_partition_node[2,9] circuit_node[2] circuit_partition_node[1,0]
 circuit_partition_node[1,1] circuit_partition_node[1,2]
 circuit_partition_node[1,3] circuit_partition_node[1,4]
 circuit_partition_node[1,5] circuit_partition_node[1,6]
 circuit_partition_node[1,7] circuit_partition_node[1,8]
 circuit_partition_node[1,9] circuit_node[1] circuit_partition_node[0,0]
 circuit_partition_node[0,1] circuit_partition_node[0,2]
 circuit_partition_node[0,3] circuit_partition_node[0,4]
 circuit_partition_node[0,5] circuit_partition_node[0,6]
 circuit_partition_node[0,7] circuit_partition_node[0,8]
 circuit_partition_node[0,9] circuit_node[0] circuit_arc[4,5]
 circuit_arc[5,4] circuit_arc[4,6] circuit_arc[6,4] circuit_arc[4,7]
 circuit_arc[7,4] circuit_arc[4,8] circuit_arc[8,4] circuit_arc[4,9]
 circuit_arc[9,4] circuit_arc[3,9] circuit_arc[9,3] circuit_arc[3,8]
 circuit_arc[8,3] circuit_arc[3,7] circuit_arc[7,3] circuit_arc[3,6]
 circuit_arc[6,3] circuit_arc[3,5] circuit_arc[5,3] circuit_arc[3,4]
 circuit_arc[4,3] circuit_arc[5,6] circuit_arc[6,5] circuit_arc[5,7]
 circuit_arc[7,5] circuit_arc[5,8] circuit_arc[8,5] circuit_arc[5,9]
 circuit_arc[9,5] circuit_arc[2,9] circuit_arc[9,2] circuit_arc[2,8]
 circuit_arc[8,2] circuit_arc[2,7] circuit_arc[7,2] circuit_arc[2,6]
 circuit_arc[6,2] circuit_arc[2,5] circuit_arc[5,2] circuit_arc[2,4]
 circuit_arc[4,2] circuit_arc[2,3] circuit_arc[3,2] circuit_arc[6,7]
 circuit_arc[7,6] circuit_arc[6,8] circuit_arc[8,6] circuit_arc[6,9]
 circuit_arc[9,6] circuit_arc[1,9] circuit_arc[9,1] circuit_arc[1,8]
 circuit_arc[8,1] circuit_arc[1,7] circuit_arc[7,1] circuit_arc[1,6]
 circuit_arc[6,1] circuit_arc[1,5] circuit_arc[5,1] circuit_arc[1,4]
 circuit_arc[4,1] circuit_arc[1,3] circuit_arc[3,1] circuit_arc[1,2]
 circuit_arc[2,1] circuit_arc[7,8] circuit_arc[8,7] circuit_arc[7,9]
 circuit_arc[9,7] circuit_arc[0,9] circuit_arc[9,0] circuit_arc[0,8]
 circuit_arc[8,0] circuit_arc[0,7] circuit_arc[7,0] circuit_arc[0,6]
 circuit_arc[6,0] circuit_arc[0,5] circuit_arc[5,0] circuit_arc[0,4]
 circuit_arc[4,0] circuit_arc[0,3] circuit_arc[3,0] circuit_arc[0,2]
 circuit_arc[2,0] circuit_arc[0,1] circuit_arc[1,0] circuit_arc[8,9]
 circuit_arc[9,8] star_edge[4,5] star_edge[4,6] star_edge[4,7]
 star_edge[4,8] star_edge[4,9] star_edge[3,9] star_edge[3,8] star_edge[3,7]
 star_edge[3,6] star_edge[3,5] star_edge[3,4] star_edge[5,6] star_edge[5,7]
 star_edge[5,8] star_edge[5,9] star_edge[2,9] star_edge[2,8] star_edge[2,7]
 star_edge[2,6] star_edge[2,5] star_edge[2,4] star_edge[2,3] star_edge[6,7]
 star_edge[6,8] star_edge[6,9] star_edge[1,9] star_edge[1,8] star_edge[1,7]
 star_edge[1,6] star_edge[1,5] star_edge[1,4] star_edge[1,3] star_edge[1,2]
 star_edge[7,8] star_edge[7,9] star_edge[0,9] star_edge[0,8] star_edge[0,7]
 star_edge[0,6] star_edge[0,5] star_edge[0,4] star_edge[0,3] star_edge[0,2]
 star_edge[0,1] star_edge[8,9]
Generals
 circuit_order[9] circuit_order[8] circuit_order[7] circuit_order[6]
 circuit_order[5] circuit_order[4] circuit_order[3] circuit_order[2]
 circuit_order[1] circuit_order[0]
End
