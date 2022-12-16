## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)

seed <- arimasim(a = 280000,  z = 281000, n = 10, p = 1, d = 0, q = 0, ar11 = 0.8, sd = 1, j1 = 4, arr1 = "0.80", n_cores = 1)
usethis::use_data(seed, internal = TRUE, overwrite = TRUE)
