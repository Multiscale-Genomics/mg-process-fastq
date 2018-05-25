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
    make_option(c("-f", "--file"), type="character", default="",
                help="Sleuth R object", metavar="character"),
    make_option(c("-t", "--tag"), type="character", default="",
                help="Tag to save the images with", metavar="character"),
    make_option(c("-l", "--degl"), type="character", default=0.05,
                help="Cutoff for differentially expressed genes", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

so = sleuth_load(file.path(opt$file))
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= opt$degl)

png(file.path(paste(opt$file, "_sample_heatmap_", opt$tag, ".png", sep="")))
plot_sample_heatmap(so)
dev.off()

png(file.path(paste(opt$file, "_transcript_heatmap_", opt$tag, ".png", sep="")), width=1000, height=2000)
plot_transcript_heatmap(so, sleuth_significant$target_id[0:100], 'est_counts')
dev.off()
