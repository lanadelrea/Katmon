final String contig = "MN908947.3";
final int[] mutpos = {670, 2445, 2790, 2832, 4184, 5406, 6770, 6786, 8208, 8991, 9344, 9534, 9866, 10029, 10449, 11537, 12030, 16616, 17410, 17678, 18163, 19955, 21618, 22020, 22054, 22200, 22295, 22332, 22578, 22599, 22674, 22679, 22686, 22688, 22775, 22786, 22813, 22882, 22893, 22898, 22910, 22942, 22992, 23012, 23012, 23013, 23013, 23013, 23013, 23018, 23040, 23055, 23063, 23075, 23525, 23599, 23604, 23854, 23948, 24424, 24469, 24503, 25810, 26060, 26270, 26529, 26530, 26577, 26709, 27382, 27383, 27384, 28271, 28311, 28881, 28882, 28883, 29510};
final char[] mutbase = {'G', 'T', 'T', 'G', 'A', 'C', 'G', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'A', 'G', 'G', 'A', 'T', 'T', 'G', 'T', 'T', 'C', 'G', 'G', 'A', 'A', 'A', 'A', 'T', 'C', 'T', 'G', 'A', 'C', 'T', 'G', 'G', 'A', 'G', 'G', 'A', 'A', 'A', 'C', 'G', 'C', 'G', 'G', 'G', 'G', 'T', 'C', 'T', 'G', 'A', 'A', 'T', 'T', 'A', 'T', 'T', 'T', 'T', 'A', 'G', 'G', 'A', 'C', 'T', 'C', 'T', 'T', 'A', 'A', 'C', 'C'};
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
