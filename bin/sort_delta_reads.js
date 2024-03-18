final String contig = "MN908947.3";
final int[] mutpos = {1048, 1148, 1191, 3037, 4181, 5184, 6402, 7124, 9053, 9891, 11201, 11418, 11514, 15243, 15451, 16466, 19220, 20320, 21618, 21846, 22227, 22335, 22916, 22916, 22917, 23604, 24410, 25352, 25469, 25855, 26104, 26767, 27638, 27739, 27752, 27874, 28461, 28881, 28916, 29402, 29427, 29742};
final char[] mutbase = {'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'G', 'C', 'T', 'T', 'A', 'T', 'T', 'T', 'G', 'T', 'T', 'T', 'A', 'G', 'A', 'G', 'A', 'T', 'T', 'T', 'T', 'C', 'C', 'T', 'T', 'T', 'G', 'T', 'T', 'T', 'A', 'T'};
int arrayLength = mutbase.length;
if(record.getReadUnmappedFlag())  return false;
if(!record.getContig().equals(contig)) return false;
ArrayList<htsjdk.samtools.SAMRecord> ret = new ArrayList<htsjdk.samtools.SAMRecord>();
for (int i=0; i < arrayLength; i++){
  if(record.getEnd() < mutpos[i]) continue;
  if(record.getStart() > mutpos[i]) continue;
  int readpos = record.getReadPositionAtReferencePosition(mutpos[i]);
  if(readpos<1) continue;
  readpos--;
  final byte[]    bases= record.getReadBases();
  if(bases[readpos]==mutbase[i]) {
    ret.add(record);
  }
  continue;
}
return ret;
