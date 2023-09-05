# Delete all lines before %%%%%%%%%% START, and save file
# sed -i "0,/%* START/d" calibrateMod3.tex
#
# Delete all lines after %%%%%%%%%% END, and save file
# sed -i '/%* END/,$d' calibrateMod3.tex

# See also: https://tex.stackexchange.com/questions/136527/section-numbering-without-numbers
#
#

# My idea:
#
# have seperate .Rnw files for each section; these .Rnw should have
# numberless section numbers, and be self contatined. That's what the
# link is for.
#
# Once those files exist, we want to be able to combine everything into
# a single document. That's why we use the sed commands to get rid of headers /
# begin and end documents. We can just save the resulting tex file
#  to the inputs/ folder, and load that later in the si.Rnw document. That way
#  we have individual sections, but also a single document.
#
#  A make file will be necessary to make this work as desired.
