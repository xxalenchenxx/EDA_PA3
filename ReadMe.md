compile
    make clean && make
run
    ./.exe <input directory> <input file.aux> <output directory>

    example
        ./PA3 toy toy.aux output
        ./PA3 ibm01 ibm01.aux output

#<output directory>: fixed, I called it as "output" in program
#if you want to visualize placement result, run following instruction
                        python print.py
 and you will get figure called "output.jpg" to see result.

#output directory will delete all files before writing new files in program