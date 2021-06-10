#!/usr/bin/env python3

import acts
import acts.examples

u = acts.UnitConstants

# cfg = acts.examples.Sequencer.Config()
# cfg.events = 10
# cfg.numThreads = 1
# s = acts.examples.Sequencer(cfg)

s = acts.examples.Sequencer(events=1000, numThreads=20)

rnd = acts.examples.RandomNumbers(seed=42)

# print("Building generic detector")
# gdc = acts.examples.GenericDetector.Config()
# detector = acts.examples.GenericDetector()
# trackingGeometry, contextDecorators = detector.finalize(gdc, None)

dd4hepCfg = acts.examples.DD4hepDetector.Config()
dd4hepCfg.xmlFileNames = ["thirdparty/OpenDataDetector/xml/OpenDataDetector.xml"]
detector = acts.examples.DD4hepDetector()
trackingGeometry, contextDecorators = detector.finalize(dd4hepCfg, None)

for cdr in contextDecorators:
    s.addContextDecorator(cdr)

field = acts.NullBField()
# field = acts.ConstantBField(acts.Vector3(0, 0, 2 * acts.UnitConstants.T))
# solenoid = acts.SolenoidBField(
#     radius = 1200*u.mm,
#     length = 6000*u.mm,
#     bMagCenter = 2*u.T,
#     nCoils = 1194
# )
# field = acts.solenoidFieldMap(
#     rlim=(0, 1200*u.mm),
#     zlim=(-5000*u.mm, 5000*u.mm),
#     nbins=(50, 50),
#     field=solenoid
# )
nav = acts.Navigator(trackingGeometry)

stepper = acts.EigenStepper(field)
# stepper = acts.AtlasStepper(field)
# stepper = acts.StraightLineStepper()
print("We're running with:", type(stepper).__name__)
prop = acts.Propagator(stepper, nav)

# algCfg = acts.examples.PropagationAlgorithm.Config(
#     prop,
#     randomNumberSvc = rnd,
#     ntests = 1000,
#     sterileLogger = False,
#     propagationStepCollection = "propagation-steps"
# )
# alg = acts.examples.PropagationAlgorithm(algCfg, acts.logging.Level.INFO)

alg = acts.examples.PropagationAlgorithm(
    propagator=prop,
    level=acts.logging.INFO,
    randomNumberSvc=rnd,
    ntests=1000,
    sterileLogger=True,
    propagationStepCollection="propagation-steps",
)

s.addAlgorithm(alg)


# output
# s.addWriter(
#     acts.examples.ObjPropagationStepsWriter(
#         collection="propagation-steps", outputDir="obj"
#     )
# )

# s.addWriter(
#     acts.examples.RootPropagationStepsWriter(
#         collection="propagation-steps", filePath="propagation_steps_python.root"
#     )
# )


s.run()
