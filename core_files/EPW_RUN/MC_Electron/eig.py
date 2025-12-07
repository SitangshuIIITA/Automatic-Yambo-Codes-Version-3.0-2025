def process_file(input_file, output_file):
    # Step 1: Read all lines from the file
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Step 2: Remove leading spaces from each line
    lines = [line.lstrip() for line in lines]

    # Step 3: Remove the last line from the file (if it's empty or doesn't contain data)
    if lines:
        lines.pop()  # Remove the last line

    # Step 4: Split the data into blocks based on empty lines
    blocks = []
    current_block = []

    for line in lines:
        if line.strip():  # If the line is not empty
            current_block.append(line.strip())  # Add line to the current block
        else:
            if current_block:
                blocks.append(current_block)  # When we hit a blank line, save the current block
                current_block = []  # Start a new block

    # Add the last block if any
    if current_block:
        blocks.append(current_block)

    # Step 5: Process each block
    processed_data = []

    # Loop through the k-values from the first block (assuming they are the same across all blocks)
    k_values = []
    for line in blocks[0]:  # First block (energies for column 2)
        k_values.append(float(line.split()[0]))  # Extract k-value from the first column of the first block

    # Loop through the k-values, and for each k, collect the energy values from all blocks
    for k_value in k_values:
        row = [f"{k_value:.4f}"]  # Start with the k-value for this row

        # For each block, find the corresponding energy for the current k-value
        for block in blocks:
            energy_found = False
            for line in block:
                values = line.split()
                if float(values[0]) == k_value:  # If k-value matches
                    energy_found = True
                    row.append(f"{float(values[1]):.4f}")  # Append the corresponding energy
                    break
            if not energy_found:
                row.append("0.0000")  # If the k-value is not found, append a placeholder

        # Add the row for this k-value
        processed_data.append("  ".join(row) + "\n")

    # Step 6: Write the processed data to the output file
    with open(output_file, 'w') as f:
        f.writelines(processed_data)

    print(f"Processed data saved to {output_file}")

# Example usage
process_file('band.dat', 'elec_disp.dat')  # Processing band.dat -> elec_disp.dat
process_file('freq.dat', 'ph_disp.dat')  # Processing freq.dat -> ph_disp.dat
