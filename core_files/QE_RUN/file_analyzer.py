import pandas as pd
class file_analyzer:
    def __init__(self):
        pass
    
    def write_csv(self, file_path, csv_file_path):
        fw = open(csv_file_path, "w")              #creating a csv file in write mode,it will create a new empty file
        with open(file_path) as f:            # opening  the given bandx file
            for line in f:  
                s = ""                        #looping through all the rows of the bandx file which doenot start with '#'
                if '#' not in line:
                    for word in line.split():       #looping through all the values in a row,word represent one value of the row
                        s = s + word + ','  
                    if len(s) == 0:
                        continue                 
                    s = s[:-1]                   #removing the last comma 
                    s = s + '\n'                #adding the value in the string
                    fw.write(s)             #write the string in the file  
                    s = ""  
        f.close()                                   #closing all the files
        fw.close()
    
    def get_distinct_x(self, table_values):

        '''in bandx file same value for x axis is repeating many times ,so to find unique x axis values we will iterate till its start repeating'''

        inti = table_values[0, 0]
    
        for i in range(1, len(table_values)):                        #iterating through all the rows
            if table_values[i,0] == inti:           #when the element of the 0th column = first element means repeatition started and we return the index
                return i

    def get_offset(self, file_path):

        file = open(file_path, "r")

        flag = "highest occupied, lowest unoccupied level (ev):"

        highest_occupied = "123"
        lowest_unoccupied = "223"
        for line in file:
            if line.strip().startswith(flag):

                data = line.strip().split(" ")
                
                for i in range(len(data)):
                    if data[i] == '': 
                        continue

                    
                    highest_occupied = lowest_unoccupied
                    lowest_unoccupied = data[i]
                
                break

        

        return [float(highest_occupied), float(lowest_unoccupied)]
    

    def get_kpoints(self, file_path):
        result = []
        print(file_path)
        with open(file_path, 'r') as file:
            data = file.readlines()

            # print(data)

            for i in range(len(data)):
                line = data[i].strip()
                if line.startswith('high-symmetry point:'):
                    values = line.split(" ")

                    val = ''
                    for word in values:
                        val = word

                    val = float(val)
                    result.append(val)

                if line == 'Bands written to file Bandx.dat':
                    break

            result
            return result
                


                    


