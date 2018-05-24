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
    make_option(c("-f", "--file"), type="character", default=NULL,
                help="Sleuth R object", metavar="character"),
    make_option(c("-t", "--tag"), type="character", default=NULL,
                help="Tag to save the images with", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

sleuth_load(file.path(opt$file))

png(file.path(opt$file + "_sample_heatmap_" + opt$save + ".png"))
plot_sample_heatmap(so)
dev.off()

png(file.path(opt$file + "_transcript_heatmap_" + opt$save + ".png"))
plot_transcript_heatmap(so, sleuth_significant$target_id[0:50])
dev.off()
