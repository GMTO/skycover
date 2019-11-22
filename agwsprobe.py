import sys
import pylab as pl


# The arguments to the script are
#    --  buffer around the probe in millimeters. The minimum distance between probes is twice this value
#    --  diameter of probe shadow
b = int(sys.argv[1])
shadow = int(sys.argv[2])

ltype = ["-","--"]

# Run this twice so the plot shows the probe outline with no buffer (solid)
# and with the buffer (dashed)
# 
for i,buffer in enumerate([0, b]):

    baffle_tube_width =  64
    baffle_tube_front = -76
    baffle_tube_back  = 409

    baffle_tube_front -= buffer
    baffle_tube_width += buffer
    
    slider_shaft_width = 66
    slider_shaft_front = 84
    slider_shaft_back  = 544

    slider_shaft_back  -= buffer
    slider_shaft_width += buffer

    slider_body_width = 330
    slider_body_front = slider_shaft_back - 1
    slider_body_back  = 1300

    slider_body_width += buffer

    pl.plot([-baffle_tube_width, baffle_tube_width, baffle_tube_width, -baffle_tube_width, -baffle_tube_width],
            [ baffle_tube_back,  baffle_tube_back,  baffle_tube_front,  baffle_tube_front,  baffle_tube_back],
            color="b", linestyle=ltype[i]
            )

    pl.plot([-slider_shaft_width, -slider_shaft_width, -slider_body_width, -slider_body_width, slider_body_width, slider_body_width, slider_shaft_width, slider_shaft_width, -slider_shaft_width],
            [slider_shaft_front,  slider_shaft_back,    slider_body_front,  slider_body_back,  slider_body_back,  slider_body_front, slider_shaft_back,  slider_shaft_front,  slider_shaft_front],
            color="r", linestyle=ltype[i]
           )
    pl.plot(0,0, "b+")

    if (i==0):
        pl.text(-50, baffle_tube_front-b-50, str(2*baffle_tube_width),color="b")
        pl.text(-50, slider_shaft_front+20, str(2*slider_shaft_width),color="r")
        pl.text(-50, slider_body_back+20,  str(2*slider_body_width),color="r")
        pl.text(baffle_tube_width+40, baffle_tube_front, str(baffle_tube_front), color="b")
        pl.text(slider_shaft_width+40, slider_shaft_front, str(slider_shaft_front)+"+s", color="r")
        pl.text(slider_body_width+40, slider_body_front, str(slider_body_front)+"+s", color="r")

pl.axis('equal')
pl.axis('off')
pl.savefig("probe.png")

prefix = "probe_"

f = open(prefix + "baffle_tube.txt", "w+")
print >>f, "{}\t{}".format(-baffle_tube_width, baffle_tube_back)
print >>f, "{}\t{}".format( baffle_tube_width, baffle_tube_back)
print >>f, "{}\t{}".format( baffle_tube_width, baffle_tube_front)
print >>f, "{}\t{}".format(-baffle_tube_width, baffle_tube_front)
f.close()

f = open(prefix + "slider_shaft.txt", "w+")
print >>f, "{}\t{}".format(-slider_shaft_width, slider_shaft_back)
print >>f, "{}\t{}".format( slider_shaft_width, slider_shaft_back)
print >>f, "{}\t{}".format( slider_shaft_width, slider_shaft_front)
print >>f, "{}\t{}".format(-slider_shaft_width, slider_shaft_front)
f.close()

f = open(prefix + "slider_body.txt", "w+")
print >>f, "{}\t{}".format(-slider_body_width, slider_body_back)
print >>f, "{}\t{}".format( slider_body_width, slider_body_back)
print >>f, "{}\t{}".format( slider_body_width, slider_body_front)
print >>f, "{}\t{}".format(-slider_body_width, slider_body_front)
f.close()

baffle_tube_front -= (shadow/2 - buffer)
baffle_tube_width += (shadow/2 - buffer)
slider_body_front  -= (shadow/2 - buffer)
slider_shaft_width += (shadow/2 - buffer)
slider_body_width += (shadow/2 - buffer)

prefix = "shadow_"

f = open(prefix + "baffle_tube.txt", "w+")
print >>f, "{}\t{}".format(-baffle_tube_width, baffle_tube_back)
print >>f, "{}\t{}".format( baffle_tube_width, baffle_tube_back)
print >>f, "{}\t{}".format( baffle_tube_width, baffle_tube_front)
print >>f, "{}\t{}".format(-baffle_tube_width, baffle_tube_front)
f.close()

f = open(prefix + "slider_shaft.txt", "w+")
print >>f, "{}\t{}".format(-slider_shaft_width, slider_shaft_back)
print >>f, "{}\t{}".format( slider_shaft_width, slider_shaft_back)
print >>f, "{}\t{}".format( slider_shaft_width, slider_shaft_front)
print >>f, "{}\t{}".format(-slider_shaft_width, slider_shaft_front)
f.close()

f = open(prefix + "slider_body.txt", "w+")
print >>f, "{}\t{}".format(-slider_body_width, slider_body_back)
print >>f, "{}\t{}".format( slider_body_width, slider_body_back)
print >>f, "{}\t{}".format( slider_body_width, slider_body_front)
print >>f, "{}\t{}".format(-slider_body_width, slider_body_front)
f.close()


