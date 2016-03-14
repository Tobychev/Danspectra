from wipfile import *

np.seterr(invalid="raise")
s6405_t5p.normalize()
mes = {}
FeImes   = FeI.measure(s6405_t5p)

vel = 0; bot = 1; con = 2; err = 3; ew  = 4; mn  = 5; var = 6; ske = 7; kur = 8
cuts = FeImes[err] < np.percentile(FeImes[err],89)

pl.subplot(3,2,1)
pl.plot(FeImes[con][cuts],FeImes[ew][cuts],'bo',alpha=0.2);
pl.title("Line bottom , FeI " + str(FeI))
pl.ylabel("Equivalent width")
pl.xlabel("Continuum intensity")

pl.subplot(3,2,2)
pl.plot(FeImes[con][cuts],FeImes[vel][cuts],'bo',alpha=0.2);
pl.title("Line bottom , FeI " + str(FeI))
pl.ylabel("Line centre")
pl.xlabel("Continuum intensity")

pl.subplot(3,2,3)
pl.plot(FeImes[con][cuts],FeImes[ske][cuts],'bo',alpha=0.2);
pl.title("Line bottom , FeI " + str(FeI))
pl.ylabel("Skewness")
pl.xlabel("Continuum intensity")

pl.subplot(3,2,4)
pl.plot(FeImes[con][cuts],FeImes[bot][cuts],'bo',alpha=0.2);
pl.title("Line bottom , FeI " + str(FeI))
pl.ylabel("Line min intesity")
pl.xlabel("Continuum intensity")

pl.subplot(3,2,5)
pl.plot(FeImes[con][cuts],FeImes[bot][cuts]/FeImes[con][cuts],'bo',alpha=0.2);
pl.title("Line bottom , FeI " + str(FeI))
pl.ylabel("Relative Line min intesity")
pl.xlabel("Continuum intensity")

pl.subplot(3,2,6)
pl.plot(FeImes[con][cuts],FeImes[var][cuts],'bo',alpha=0.2);
pl.title("Line bottom , FeI " + str(FeI))
pl.ylabel("variance")
pl.xlabel("Continuum intensity")

pl.show()

