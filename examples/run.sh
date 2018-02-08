
## this will set the path for the runs here. Effectively, it changes
## the system variable path to include also the ../ directory, where
## the comus executable is located. Modify it if you have put the
## comus executable somewhere else
COMUS_PATH=../

## this will set the path. It has an effect only on this script
PATH=$COMUS_PATH:$PATH
echo $PATH


# Example 1

# comus as Hudson's ms. Comus can function as Hudson's ms. The output
# will be sent to ms.out. Note that simulations are performed under
# the infinite site model. For more details please read the manual of
# Hudson's ms
comus 10 10 -t 6 -r 100 1000 > ms.out


# Example 2

##-name

## when you simulate more than one species (or you specify explicitly
## 1 species) then you should use the -name flag. Thus, output files
## will be generated that will contain the results

## Infinite site model - more than 1 species
comus 2 5 6 15 -t 10 -r 10 -rsites 100 -name run1


# Example 3

# Finite site model - single species.  Simulation of a single species
# under the finite site model
comus 1 10 10 -t 6 -mm hky -msites 100 -name run2

# -r <rec rate> -rsites <number of segments that can be recombined>
# two species and recombination rate
comus 2 5 8 1 -t 20 -r 10 -rsites 100  -mm jc -msites 200 -name run3 -noseparator


# Example 4

# -b <birth rate>
# specify the birth rate for the generation of the guide tree
comus 2 5 8 1 -t 20 -birth 1 -r 10 -rsites 100  -mm jc -msites 200 -name run4 -noseparator

# You can compare results from run3 to results from run4
# e.g. run seaview comus_Results.run3 
# and then seaview comus_Results.run4

# you will see that in comus_Results.run4 the first 5 sequences are
# quit separated from the next 8. The reason is that the birth rate
# is much smaller in run4, therefore the branch length is considerably
# larger. On the other hand, the sequences in run3 cannot be separated so well. 


# Example 5

# -phylomode < 0|1|2|3 > -torigin <double> 
# The default value is 3, i.e. simulations by fixing the number of taxa

# specify either the time that the birth-death process starts or the
# time of the most recent common ancestor. This depends on the
# phylomode. 
# IMPORTANT: The flag -phylomode should be before the -torigin. 
# -torigin works only with -phylomode 1 or -phylomode 2
# Thus, with -phylomode 1 you set the TMRCA
comus 3  2 3 3  1 -t 10 -mm F84 -phylomode 1 -torigin 8.4  -oPhylo -name run5 -noSeparator

## Check the comus_Phylo.run5 to confirm that the torigin has been set
## to 8.4


# Example 6 -phylomode 2 with -phylomode 2 you set the time where the
# speciation process starts. Note that trees can be dramatically
# different if you set the -phylomode 1 or -phylomode 2. Eventually
# the shape of the tree will depend on the combination of the birth
# (default 100), death (default 0), sampling (default 1), the number
# of species and the -phylomode option. This example is the same as
# run5 but here we use -phylomode 2
comus 3  2 3 3  1 -t 10 -mm F84 -phylomode 2 -torigin 8.4  -oPhylo -name run6 -noSeparator

# if you compare the two phylogenetic guide trees, from run5 and run6
#pavlos@mira:~/research/COMUS/comus/examples$ more comus_Phylo.run6
#((1:0.0113582768669188:0, 3:0.0113582768669188:0):0.0024958592294403:2:[1 3 ], 2:0.0138541360963590:0);

#pavlos@mira:~/research/COMUS/comus/examples$ more comus_Phylo.run5
#(1:8.4000000000000004:0, (3:0.0032581650891877:0, 2:0.0032581650891877:0):8.3967418349108129:2:[3 2 ]);

# you see that in the second tree, the time of the root is 8.4,
# whereas in the first one, the time of the root is about 0.0138. In
# the second tree, the process has started at time 8.4, but the first
# speciation event (i.e. the root) happened at time 0.0138. You may
# wonder why the time of the root is so much more recent than the time
# that the process starts (8.4 vs 0.0138). To understand this, you
# should take into account the birth rate as well. By default birth
# rate is 100, i.e. quite large. If the first speciation event will
# take place soon after the process starts then, the probability that
# only one more speciation will happen (total number of species is 3)
# will happen while two branches exist is quite small. Thus, the
# scenario where no speciation happens (and therefore only one branch
# is present) till very recently is much more probable. See next
# example, the situation will be much more different when the
# speciation rate is small

#Example 7
# the same as the previous example, but the speciation rate is now 0.01
comus 3  2 3 3  1 -t 10 -mm F84 -phylomode 2 -torigin 8.4  -birth 0.01 -oPhylo -name run7 -noSeparator

# now the first speciation event happend at time  7.6
#pavlos@mira:~/research/COMUS/comus/examples$ more comus_Phylo.run7
#((1:0.2900347629188417:0, 3:0.2900347629188417:0):7.3566725603605523:2:[1 3 ], 2:7.6467073232793936:0);


#Example 8 
#-phylomode 0 
#this is the speciation process ala
#yang-rannala. It is the same as -phylomode 1 but now the TMRCA has
#been set to 1.0
comus 3  2 3 3  1 -t 10 -mm F84 -phylomode 0 -birth 0.01 -oPhylo -name run8 -noSeparator

# again when you simulate under -phylomode 0, pay attention to birth
# rate etc. Compare results run8 and run9, where the birth rate is
# small and large, respectively


# Example 9
comus 3  2 3 3  1 -t 10 -mm F84 -phylomode 0 -birth 1000 -oPhylo -name run9 -noSeparator

#::::::::::::::
#comus_Phylo.run9
#::::::::::::::
#(3:1.0000000000000000:0, (1:0.0002706493376888:0, 2:0.0002706493376888:0):0.9997293506623112:2:[1 2 ]);
#::::::::::::::
#comus_Phylo.run8
#::::::::::::::
#(1:1.0000000000000000:0, (3:0.2018062378403698:0, 2:0.2018062378403698:0):0.7981937621596302:2:[3 2 ]);

# You see that in both cases the TMRCA is at 1.0. Also, the next
# speciation event, occurred very after a long time (ie. recently,
# backwards in time) when the speciation rate is large (run9)

# Example 10
#-phylomode 4 -oldestOrigin <double>

# By using -phylomode 4 you can define the oldest time that the
# process starts.  Note that this is done by a rejection method. If
# however, the algorithm tries to simulate such a genealogy without
# success for a large number of trials (>10000), then it throws an
# error For example, the following command where the birth rate is
# small and the oldest origin very young will probably result in an
# error:
# comus: comus_phylo.c:1086: BranchLengthBD_Stadler_ntaxa_ProcessOldestStart: Assertion `cnt < 10000' failed.
# ./run.sh: line 156:  7394 Aborted 
comus 3  2 3 3  1 -t 10 -mm F84 -phylomode 4 -oldestOrigin 0.01 -birth 0.001 -oPhylo -name run10 -noSeparator

# Example 11
#In contrast, this command, where the speciation rate is larger, will probably run just fine
comus 3  2 3 3  1 -t 10 -mm F84 -phylomode 4 -oldestOrigin 0.01 -birth 15 -oPhylo -name run11 -noSeparator


#Example 12,13,14
#-death This flag specifies the death rate and of course it
#affects the phylogenetic tree simulation. The death rate can be
#greater than the birth rate. In this case however, when death rate is
#greater than birth rate, it does NOT affect the process of
#speciation. Please note that the differences in birth and death
#rates, do not affect dramatically the process. This is probably
#because the times of the nodes are given by:

# t[i] = 1/(lamb1-mu1)*log((lamb1-mu1* exp((-lamb1+mu1)*torigin) -mu1*(1-exp((-lamb1+mu1)*torigin)) *r )/(lamb1\
#          -mu1* exp((-lamb1+mu1)*torigin) - lamb1*(1-exp((-lamb1+mu1)*torigin)) *r )   )  ;

# thus, time depends mostly on the product (lamb1 - mu1 ) * torigin. 
# when the birth - death ~ 0, then torigin is large and lamb1 - mu1 is small. On the other hand, when birth is large compared to death, then torigin is small, and lamb1 - mu1 is large. Check examples 12 and 13, 14. 

comus 3 2 3 3 1 -t 10 -mm F84  -birth 15000 -death 14999.999 -oPhylo -name run12 -noSeparator

comus 3 2 3 3 1 -t 10 -mm F84  -birth 15000 -death 0 -oPhylo -name run13 -noSeparator

comus 3 2 3 3 1 -t 10 -mm F84  -birth 15000 -death 14999.999  -phylomode 2 -torigin 1 -oPhylo -name run14 -noSeparator


# Specifying demographic changes within species boundaries

# Example 15 

#-eN <time> <species> <newsize>

# with this command you can change the size of all populations of one species

# As in ms, we can specify demographic changes within the species boundaries. The next model, specifies a population decline (forward in time) for species 1, or population expansion backward in time. Obviously, since species 1 was larger, the amount of variability should be larger within species 1. You can inspect that by either viewing the comus_Results.run15, or the file with the coalescent trees comus_Coalescent.run15 (see that the branch lengths associated with species 1 are larger). Also you can build a phylogenetic tree using raxml and view the inferred branch lengths. For species 1 the branches will be larger

comus 4 10 10 10 10 1 -t 10 -msites 10000 -mm JC -birth 0.1 -eN 0.000001 1 1000 -oPhylo -oCoalescent -noSeparator -name run15

#with the next two commands you can build the phylogenetic tree using
#raxml (you need to install it first). Then you can view the tree
#using njplot

# raxmlHPC -n raxTree15 -s comus_Results.run15 -m GTRGAMMA
# njplot RAxML_bestTree.raxTree15 




# Example 16

#-eN <time> <new size>

#with this command you can change the sizes of all populations simultaneously

comus 4 10 10 10 10 1 -t 10 -msites 10000 -mm JC -birth 0.1 -eN 0.000001 1000 -oPhylo -oCoalescent -noSeparator -name run16


# Example 17, 18
#-noSeparator (case insensitive)

#for many of the previous examples, we have used the
#-noSeparator. With -noSeparator the // flag is not printed. Thus, we
#can open the alignment result easier in an alignment viewer. For
#example, compare the comus_Results files of run17 and run18
comus 4 10 10 10 10 1 -t 10 -msites 10000 -mm JC -birth 0.1 -eN 0.000001 1000 -oPhylo -oCoalescent -noSeparator -name run17

comus 4 10 10 10 10 1 -t 10 -msites 10000 -mm JC -birth 0.1 -eN 0.000001 1000 -oPhylo -oCoalescent -name run18


# Example 19 gradual isolation model.  with the gradual isolation
# model, we allow gene flow to occur between two recently generated
# species. The file tree.newick2 specifies that speciation occurred at
# time 1.0. However, the species continue to be in partial contact for
# a period 0.8. Thus, it is possible that coalescent between lineages
# of the two species will occur more recently than speciation. And
# this is what probably will happen in run20 (of course the process is
# stochastic, so it is possible that it will not happen). However, the
# probability should be large. So if you open the
# comus_Coalescent.run20 you should see several coalescent trees with
# coalescent more recent than the speciation. Units of partIsoPeriod
# are phylogenetics
comus 2 1 1 1 -iphylo tree.newick2 -t 10 -mm hky -name run19 -oPhylo -oCoalescent
comus 2 1 1 100 -iphylo tree.newick2 -t 10 -mm hky -name run20 -oPhylo -oCoalescent -partIsolation -partMaxMigration 200 -partIsoPeriod 0.8



#Ancestral sampling
#comus can do ancestral sampling on a population-based level.
# this means that it is possible to sample the sequences of some population (of some species) at some point in the past.
# This is useful when the study includes fossils, for example a common study between neandertals and h. sapiens. 
# How this is implemented:
#Let's assume that sequences from population A have been sampled at
#time point t0 (in phylogenetic units) in the past. The usual time
#conversion takes place,  
#speciation events are translated to population merge events (as
#usually). However now, all events that include population A are NOT
#allowed till moment t0. Effectively this generates trees that look
#like:

#                                                                        +---+ 9 
#                                                                        |       
# +----------------------------------------------------------------------+ +-+ 6 
# |                                                                      | |     
# |                                                                      +-+++ 10
# |                                                                        ||    
# |                                                                        +++ 7 
# |                                                                         |    
#=|                                                                         ++ 8 
# |                                                                              
# |                                    ++ 3                                      
# |                                    |                                         
# |                                   ++| 1                                      
# |                                   |++                                        
# +-----------------------------------+ | 4                                      
#                                     |                                          
#                                     |-+ 2                                      
#                                     |                                          
#                                     +-+ 5                                      
                                                              

# The coalescent tree above has been generated with the following command

# Example 21 
comus 2 5 5 1   -t 1000 -mm hky -msites 10000  -name run21 -oPhylo -oCoalescent -iphylo tree.newick2 -samplingTime 1 0.5


## Rate Heterogenity
comus 2 5 5 1 -t 100 -mm hky -msites 10000 -name run22 -gammaCategories 4 -alpha 0.1 -iphylo tree.newick2 -samplingTime 1 0.5
