// Read mutation file path from system property
String mutationFilePath = System.getProperty("mutationFilePath");
if (mutationFilePath == null || mutationFilePath.isEmpty()) {
    throw new IllegalArgumentException("mutationFilePath system property must be provided.");
}

// Read mutations from the file
java.util.List<Integer> positions = new java.util.ArrayList<>();
java.util.List<Character> bases = new java.util.ArrayList<>();

try (java.io.BufferedReader reader = new java.io.BufferedReader(new java.io.FileReader(mutationFilePath))) {
    String line;
    while ((line = reader.readLine()) != null) {
        line = line.trim();
        if (line.isEmpty()) continue; // Skip empty lines

        String[] parts = line.split("\\s+");
        if (parts.length != 2) {
            throw new IllegalArgumentException("Invalid line format in mutations file: " + line);
        }

        positions.add(Integer.parseInt(parts[0]));
        bases.add(parts[1].charAt(0));
    }
} catch (java.io.IOException e) {
    throw new RuntimeException("Error reading mutation file: " + e.getMessage(), e);
}

// Convert lists to arrays
int[] mutpos = positions.stream().mapToInt(i -> i).toArray();
char[] mutbase = new char[bases.size()];
for (int i = 0; i < bases.size(); i++) {
    mutbase[i] = bases.get(i);
}

// Set the contig (will modify this to be dynamic if needed)
String contig = "MN908947.3";

// Early return if the record is unmapped or doesn't match the contig
if (record.getReadUnmappedFlag() || !record.getContig().equals(contig)) {
    return null; // Skip this record
}

// Check for mutations
for (int i = 0; i < mutpos.length; i++) {
    if (record.getEnd() < mutpos[i] || record.getStart() > mutpos[i]) {
        continue; // Skip if the mutation position is outside the record's range
    }

    int readpos = record.getReadPositionAtReferencePosition(mutpos[i]);
    if (readpos < 1) {
        continue; // Skip if the position is invalid
    }

    // Check if the base at the read position matches the mutation base
    byte[] basesArray = record.getReadBases();
    if (basesArray[readpos - 1] == mutbase[i]) {
        return record; // Return the record if it matches the mutation
    }
}

return null; // Skip this record