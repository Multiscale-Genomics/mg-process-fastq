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
    make_option(c("-f", "--files"), type="character", default=NULL,
                help="Comma separated sample files", metavar="character"),
    make_option(c("-bg", "--bg_files"), type="character", default=NULL,
                help="Comma separated background files", metavar="character"),
    make_option(c("-s", "--species"), type="character", default=NULL,
                help="Species", metavar="character"),
    make_option(c("-a", "--assembly"), type="character", default=NULL,
                help="Assembly", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="BigWig file", metavar="character"),
    make_option(c("-l", "--local_lib"), type="character", default=NULL,
                help="Local R library", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

library(
    paste('BSgenome.', opt$species, '.', opt$assembly, sep=""),
    lib.loc=opt$local_lib,
    character.only=TRUE)
genome <- get(paste('BSgenome.', opt$species, '.', opt$assembly, sep=""))

fragments <- getDpnIFragments(genome)

sample_files = unlist(strsplit(opt$files, ","))
bg_files = unlist(strsplit(opt$bg_files, ","))

samples <- data.frame(name = factor(c(rep(opt$sample_name, length(sample_files)),
                                      rep(opt$background_name, length(bg_files))),
                                    levels = c(opt$sample_name, opt$background_name)),
                      replicate = c(seq(1,length(sample_files)), seq(1,length(bg_files))),
                      paths = c(sample_files, bg_files))

results <- getEnrichedDamRegions(samples,fragments)

poi.positive <- results[mcols(results)$log2FoldChange > 0, ]
poi.negative <- results[mcols(results)$log2FoldChange < 0, ]

saveBigWigScore(samples,genome,opt$output)
