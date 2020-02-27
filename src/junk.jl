using Printf
using Statistics

s = 0
s = "dajo"

println(s)

c1 = 'c'

println(c1)

i1 = UInt8(trunc(2.11))
print(i1)

# comentar

long_string = """ this
ahoj
dlouhy
long_string
jep"""


println(long_string)

if 12 > 1
    print("is bigger.")
else
    println("not biiger")

end

println()

@printf("print string %s\n", "THIS IS A STRING")

ar = rand(0:9, 3, 3)

println(ar)

function one(ahoj = 1)
    println(ahoj)
end


one(12)

one()




 