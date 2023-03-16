library(slimr)

slim_script(
  slim_block(initialize(),
             {
               initializeTreeSeq();
         initializeMutationRate(0);  #no non-selected mutations
         initializeMutationType("m1", 0.5, "f", 0.0); #doesn't matter, never will be used.
         initializeMutationType("m2", 0.5, "f", 0.01);  # POSITION1 introduced mutation, adjust s and dominance
         m2.convertToSubstitution = F;
         initializeGenomicElementType("g1", m1, 1.0);
         initializeGenomicElement(g1, 0, 20);
         initializeRecombinationRate(0);
}),
slim_block("s3",1,early(), { sim.addSubpop("p1", 15000);}), # //POSITION2 adjust population size, this appears to be diploid size
slim_block("s7", 3, late() ,{ 
         target = sample(p1.genomes, 1);
         target.addNewDrawnMutation(m2, 10);
}),
slim_block("s1",5, 1000000, early(), {
	mut = sim.mutationsOfType(m2);
	if (size(mut) == 1) { #//if the mutation is still segregating
     if  (sim.mutationFrequencies(NULL, mut) > 0.745  & sim.mutationFrequencies(NULL, mut) < 0.755) #//POSITION3
            {
                cat(sim.mutationFrequencies(NULL, mut));
              cat(": ESTABLISHED\n");  sim.treeSeqOutput("CLUES2/CLUES2SimulationStudy/TrueTrees/slim.trees"); 
              }

     if  (sim.mutationFrequencies(NULL, mut) == 1.0)
            {  cat("1.0");
               sim.simulationFinished();
              } }
     	if (size(mut) == 0)  {  cat("0.0");
               sim.simulationFinished();
              }
}),
slim_block("s99", 1000000,early(), {
	mut = sim.mutationsOfType(m2);
cat(sim.mutationFrequencies(NULL, mut) );
sim.simulationFinished();
})) -> script1

slim_run(script1, show_output = T, keep_all_output = T, simple_run = T)