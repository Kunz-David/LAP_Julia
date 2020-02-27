using Printf

@printf("Ahoj")


using PkgTemplates
t = Template(user="David-Kunz",
             authors=["David Kunz"],
             julia_version=v"1.3.1",
             dir="/Users/MrTrololord/Google_Drive/cvut/bakalarka",
             ssh=true
             )

generate(t, "LAP_julia")
