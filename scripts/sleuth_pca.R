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
    make_option(c("-v", "--covariant"), type="character", default="",
                help="Specific covariant to plot images for", metavar="character"),
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

so = sleuth_load(file.path(opt$file))

png(file.path(opt$file + "_pca_" + opt$tag + ".png"))
plot_pca(so, color_by = opt$covariant)
dev.off()

png(file.path(opt$file + "_pca1_" + opt$tag + ".png"))
plot_loadings(so, pc_input = 1)
dev.off()
