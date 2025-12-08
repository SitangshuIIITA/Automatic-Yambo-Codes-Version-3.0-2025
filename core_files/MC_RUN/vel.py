def process_file(input_file):
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

    # Step 5: Process the electron data from the relevant block
    processed_data = {}

    # Iterate through the blocks to extract the electron velocity data
    for block in blocks:
        for line in block:
            # Skip header lines or malformed lines
            if line.startswith("#") or line.strip() == "":
                continue
            
            values = line.split()
            
            # Skip lines that don't have enough columns
            if len(values) < 7:
                print(f"Skipping malformed line: {line}")
                continue
            
            # Extract k-point, band index, velocity components, eigenvalue, and weight
            try:
                ik = int(values[0])  # k-point index
                ibnd = int(values[1])  # Band index
                velocity_x = float(values[2])  # x-component of velocity
                velocity_y = float(values[3])  # y-component of velocity
                velocity_z = float(values[4])  # z-component of velocity
                eig = float(values[5])  # Eigenvalue or other associated property
                weight = float(values[6])  # Weight associated with the state

                # Store the parsed data in a dictionary with band index as key
                if ibnd not in processed_data:
                    processed_data[ibnd] = []

                processed_data[ibnd].append([ik, ibnd, velocity_x, velocity_y, velocity_z, eig, weight])

            except ValueError as e:
                print(f"Error processing line: {line}. Error: {e}")
                continue  # Skip malformed lines

    # Step 6: Write the processed data to separate files for each band
    for ibnd, data in processed_data.items():
        output_file = f"band_{ibnd}.txt"
        with open(output_file, 'w') as f:
            for row in data:
                row_str = [f"{row[0]}", f"{row[1]}", 
                           f"{row[2]:.4e}", f"{row[3]:.4e}", 
                           f"{row[4]:.4e}", f"{row[5]:.4e}", 
                           f"{row[6]:.4e}"]
                f.write("  ".join(row_str) + "\n")
        print(f"Processed data for band {ibnd} saved to {output_file}")

# Example usage
process_file('electron_IBTEvel_sup_0.fmt')

