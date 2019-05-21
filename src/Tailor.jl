global grades
global step
grades = 0
step = 0


input = ARGS[1]
output = ARGS[2]
countmax = Meta.parse(ARGS[3])

open(input) do inputfile
    lines = readlines(inputfile)

    open(output) do outputfile
        for ln in lines
            if ln[1:4] == "grades"
                eval(Meta.parse(ln))
            end

            if ln[1:4] == "step"
                eval(Meta.parse(ln))
                stepflag = 0
                recordflag = 0
            end

            if (step-1)%countmax == 0
                stepflag = 1
            end

            if ln[1:4] == "Real"
                recordflag = 1
            end

            if stepflag == 1 && recordflag ==1
                write(outputfile,"$ln\n")
            end
        end
    end

end

