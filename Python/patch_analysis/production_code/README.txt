Patch Analysis Pipeline v1.0
12.14.2012
by Mark Evans

Step 1: Generate patch files for 3D model to be used
genSurfPatches.py
- should only ever need to run once. Patches can then be reused indefinitely

Step 2: Prepare patch/IC50 vector input files
batch_prep7.py which calls pepa_seq_prep7.5.py

Step 3: Process the patch vector files through R models
processPatches_v10.py

Step 4: Create index of result files so you know what is what
create_index.py

Step 5: Create linear alignments to determine putative epitope regions
linear_align_batch.py

Step 6: Determine final epitope region and generate images and epitope lists
calcFinalEpitope2.py which calls chimeraFinalImageGenerator.py

Step 7: Optional - If you are scanning a variety of conditions, this will make HTML summaries of results
makeSummaries2.py
