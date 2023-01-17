for a in {6..15};
do
for b in {1..10};
do
MODEL=$a
NAME=_replicate$b
mv ~/GhostPopulationSimulation/Models/Model$MODEL/replicate$b/model$MODEL$NAME.vcf.2pop.pi.windowed.pi ~/GhostPopulationSimulation/Models/Model$a/PI/;
mv ~/GhostPopulationSimulation/Models/Model$MODEL/replicate$b/model$MODEL$NAME.vcf.3pop.pi.windowed.pi ~/GhostPopulationSimulation/Models/Model$a/PI/;
mv ~/GhostPopulationSimulation/Models/Model$MODEL/replicate$b/model$MODEL$NAME.vcf.2pop.fst.windowed.weir.fst ~/GhostPopulationSimulation/Models/Model$a/FST;
mv ~/GhostPopulationSimulation/Models/Model$MODEL/replicate$b/model$MODEL$NAME.vcf.3pop.fst.windowed.weir.fst ~/GhostPopulationSimulation/Models/Model$a/FST; 
done;
done