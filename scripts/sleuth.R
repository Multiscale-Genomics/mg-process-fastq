#!/usr/bin/env Rscript

# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

library("optparse")
library("sleuth")

option_list = list(
    make_option(c("-c", "--config"), type="character", default=NULL,
                help="Sleuth configuration file", metavar="character"),
    make_option(c("-o", "--save"), type="character", default=NULL,
                help="Location to save the sleuth object to", metavar="character"),
    make_option(c("-d", "--data_dir"), type="character", default=NULL,
                help="Data directory for the samples and kallisto results", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

sample_id <- dir(file.path(opt$data_dir))
kal_dirs <- file.path(opt$data_dir, sample_id, "kallisto")

s2c <- read.table(opt$config, header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample, condition)
s2c <- dplyr::mutate(s2c, path = kal_dirs)

so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, num_cores = 4)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

models(so)

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)

sleuth_save(so, opt$save)
