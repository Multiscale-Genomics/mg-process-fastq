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
library(BSgenome)
library(iDEAR)

option_list = list(
    make_option(c("-n", "--sample_name"), type="character", default=NULL,
                help="Sample Name", metavar="character"),
    make_option(c("-bn", "--background_name"), type="character", default=NULL,
                help="Background Sample Name", metavar="character"),
    make_option(c("-f1", "--file1"), type="character", default=NULL,
                help="Sample 1", metavar="character"),
    make_option(c("-f2", "--file2"), type="character", default=NULL,
                help="Sample 2", metavar="character"),
    make_option(c("-f3", "--file3"), type="character", default=NULL,
                help="Background Sample 1", metavar="character"),
    make_option(c("-f4", "--file4"), type="character", default=NULL,
                help="Background Sample 2", metavar="character"),
    make_option(c("-s", "--species"), type="character", default=NULL,
                help="Species", metavar="character"),
    make_option(c("-a", "--assembly"), type="character", default=NULL,
                help="Assembly", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="BigWig file", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

library(paste('BSgenome.', opt$species, '.', opt$a, sep=""), character.only=TRUE)
genome <- get(paste('BSgenome.', sp, '.', a, sep=""))

fragments <- getDpnIFragments(genome)

samples <- data.frame(name = factor(c(opt$sample_name, opt$sample_name, opt$background_name, opt$background_name),
                                    levels = c(opt$sample_name, opt$background_name)),
                      replicate = c(1,2,1,2),
                      paths = c(opt$file1, opt$file2, opt$file3, opt$file4))

results <- getEnrichedDamRegions(samples,fragments)

poi.positive <- results[mcols(results)$log2FoldChange > 0, ]
poi.negative <- results[mcols(results)$log2FoldChange < 0, ]

saveBigWigScore(samples,genome,opt$output)
