# Author:    Lerato E. Magosi
# R version: 4.1.2 (2021-11-01)
# Platform:  x86_64-apple-darwin17.0 (64-bit)
# Date:      16May2022
# Purpose:   eLife protocol request. Permutation test to evaluate whether there was preferential sexual mixing among trial communities by geographic proximity
# License:   MIT



# Goal: Perform a permutation test to assess whether there is preferential genetic linkage 
#       by geographic distance


# Load libraries
# -----------------------

library(dplyr)
library(tidyr)


# Load data
# -----------------------

# Load dataset of all possible bcpp community pairs excluding same community pairs and community pairs where genetic links were observed
df_dist_time_excl_obs_pairs <- read.csv("eLife_bio_protocol_request_permutation_testing_datsets/drive_distance_time_exclude_same_community_and_observed_transmission_pairs.csv")

# View dataset structure
df_dist_time_excl_obs_pairs %>% glimpse()


# Load dataset of bcpp community pairs where genetic links were observed
df_dist_time_obs_trm_pairs_with_replicates <- read.csv("eLife_bio_protocol_request_permutation_testing_datsets/drive_distance_time_observed_transmission_pairs_with_replicates.csv")

# View dataset structure
df_dist_time_obs_trm_pairs_with_replicates %>% glimpse()

# View dataset
df_dist_time_obs_trm_pairs_with_replicates


# Extract names of all 30 bcpp communities
bcpp_comm_pre_shuffle <- df_dist_time_excl_obs_pairs %>% pull(origin) %>% unique()
bcpp_comm_pre_shuffle

# Get corresponding indices for origin and destination communities in: df_dist_time_obs_trm_pairs_with_replicates from: bcpp_comm_pre_shuffle 
idx_obs_origin_comms <- lapply(df_dist_time_obs_trm_pairs_with_replicates %>% pull(origin), function(x) { output <- which(bcpp_comm_pre_shuffle == x); output })

idx_obs_origin_comms <- unlist(idx_obs_origin_comms)

idx_obs_origin_comms


# destination communities
idx_obs_destination_comms <- lapply(df_dist_time_obs_trm_pairs_with_replicates %>% pull(destination), function(x) { output <- which(bcpp_comm_pre_shuffle == x); output })

idx_obs_destination_comms <- unlist(idx_obs_destination_comms)

idx_obs_destination_comms


# Sanity check! Are we getting the right indices?
# We expect the indices obtained in the vectors: idx_obs_origin_comms ans idx_obs_destination_comms to 
# correspond the the positions of bcpp communities in the vector: bcpp_comm_pre_shuffle as named in 
# the origin and destination fields of the data.frame: df_dist_time_obs_trm_pairs_with_replicates
#
# Example: 
# The first 3 communities in the origin and destination fields of df_dist_time_obs_trm_pairs_with_replicates are:

# > df_dist_time_obs_trm_pairs_with_replicates %>% select(origin, destination)
#          origin   destination
# 1        Gumare       Shakawe
# 2  Lentsweletau     Mmankgodi
# 3       Masunga        Sebina

df_dist_time_obs_trm_pairs_with_replicates %>% select(origin, destination)

# The origin communities are located at positions: 3, 5, 8 and the destination communities at: 27, 14, 24 

# > bcpp_comm_pre_shuffle
#  [1] "Bokaa"          "Digawana"       "Gumare"         "Gweta"         
#  [5] "Lentsweletau"   "Lerala"         "Letlhakeng"     "Masunga"       
#  [9] "Mathangwane"    "Maunatlala"     "Metsimotlhabe"  "Mmadinare"     
# [13] "Mmandunyane"    "Mmankgodi"      "Mmathethe"      "Molapowabojang"
# [17] "Nata"           "Nkange"         "Oodi"           "Otse"          
# [21] "Rakops"         "Ramokgonami"    "Ranaka"         "Sebina"        
# [25] "Sefhare"        "Sefophe"        "Shakawe"        "Shoshong"      
# [29] "Tati_siding"    "Tsetsebjwe"    

bcpp_comm_pre_shuffle

# Thus we expect the first three entries of idx_obs_origin_comms and idx_obs_destination_comms to be:
# 3, 5, 8 and 27,14 and 24 respectively

idx_obs_origin_comms[1:3]

idx_obs_destination_comms[1:3]

# End of exampple


# Permute the order of communities in: bcpp_comm_pre_shuffle over a set number of iterations

# Get number of communities
number_bcpp_comms <- length(bcpp_comm_pre_shuffle)
number_bcpp_comms

# permute ordering of communities
iterations <- 10000
lst_iterations_perm_bcpp_comm <- lapply(seq(iterations), function(x) {print(paste0("iteration: ", x)); output <- bcpp_comm_pre_shuffle %>% sample(number_bcpp_comms); output } )

lst_iterations_perm_bcpp_comm


# Assemble permuted origin destination community pairings

# Get permuted origins 
# srtategy: for each iteration i.e. dataset in: lst_iterations_perm_bcpp_comm get corresponding communities for indicies provided in: idx_obs_origin_comms
perm_origins <- mapply( function(vec_iteration_perm_comms) { unlist(lapply(idx_obs_origin_comms, function(x) { vec_iteration_perm_comms[x] } )) }, lst_iterations_perm_bcpp_comm, SIMPLIFY = FALSE)
perm_origins

# Get permuted destinations
perm_destinations <- mapply( function(vec_iteration_perm_comms) { unlist(lapply(idx_obs_destination_comms, function(x) { vec_iteration_perm_comms[x] } )) }, lst_iterations_perm_bcpp_comm, SIMPLIFY = FALSE)
perm_destinations


# Combine into list of data.frames
lst_perm_orig_dest_comms <- lapply(seq_along(lst_iterations_perm_bcpp_comm), function(x) { output <- data.frame("iteration_no" = x, "origin" = perm_origins[[x]], "destination" = perm_destinations[[x]]); output <- inner_join(output, df_dist_time_excl_obs_pairs); output })

names(lst_perm_orig_dest_comms) <- paste0("iteration_", seq(iterations))

lst_perm_orig_dest_comms %>% head()

lst_perm_orig_dest_comms %>% tail()


# Combine data.frames into a single data.frame
df_combine_perm_orig_dest_comms <- dplyr::bind_rows(lst_perm_orig_dest_comms, .id = "id")

df_combine_perm_orig_dest_comms %>% glimpse()


# Get mean, min and max distance and travel time for observed between community trm pairs
# ------------------

df_mean_min_max_dist_time_obs_trm_pairs <- df_dist_time_obs_trm_pairs_with_replicates %>%
    dplyr::summarise(mean_dist_km = mean(curr_travel_dist_km),
                     min_dist_km = min(curr_travel_dist_km),
                     max_dist_km = max(curr_travel_dist_km),
                     mean_time_h = mean(curr_travel_time_h),
                     min_time_h = min(curr_travel_time_h),
                     max_time_h = max(curr_travel_time_h))

df_mean_min_max_dist_time_obs_trm_pairs


# Compare mean travel distances and times for permuted and observed community pairs

ELSE <- TRUE

df_mean_perm_orig_dest_comms <- df_combine_perm_orig_dest_comms %>% 
    dplyr::group_by(id) %>% 
    dplyr::summarise("mean_travel_dist_km_perm" = mean(curr_travel_dist_km), "mean_travel_time_h_perm" = mean(curr_travel_time_h)) %>% 
    dplyr::mutate(
           mean_travel_dist_km_obs = df_mean_min_max_dist_time_obs_trm_pairs %>% pull(mean_dist_km), 
           mean_travel_time_h_obs = df_mean_min_max_dist_time_obs_trm_pairs %>% pull(mean_time_h), 
           mean_travel_dist_km_perm_lt_obs = case_when(mean_travel_dist_km_perm < mean_travel_dist_km_obs ~ 1, ELSE ~ 0), 
           mean_travel_time_h_perm_lt_obs = case_when(mean_travel_time_h_perm < mean_travel_time_h_obs ~ 1, ELSE ~ 0)
           )

df_mean_perm_orig_dest_comms


# Count number of iterations where mean travel distances and times from permuted community pairs 
#     were less than those from observed community pairs
df_mean_perm_orig_dest_comms %>% count(mean_travel_dist_km_perm_lt_obs)

df_mean_perm_orig_dest_comms %>% count(mean_travel_time_h_perm_lt_obs)


# Compute p-values

# travel time
number_iterations_mean_traveltime_h_perm_lt_obs <- df_mean_perm_orig_dest_comms %>% 
    count(mean_travel_time_h_perm_lt_obs) %>% 
    filter(mean_travel_time_h_perm_lt_obs == 1) %>% 
    pull(n)
    
pvalue_1_sided_mean_travel_time <- number_iterations_mean_traveltime_h_perm_lt_obs / iterations
pvalue_1_sided_mean_travel_time

# travel distance
number_iterations_mean_travel_dist_km_perm_lt_obs <- df_mean_perm_orig_dest_comms %>% 
    count(mean_travel_dist_km_perm_lt_obs) %>% 
    filter(mean_travel_dist_km_perm_lt_obs == 1) %>% 
    pull(n)
    
pvalue_1_sided_mean_travel_dist <- number_iterations_mean_travel_dist_km_perm_lt_obs / iterations
pvalue_1_sided_mean_travel_dist

      

# Display session information for record keeping
# ----------------------------------------------
sessionInfo()


# And, thats all folks!


