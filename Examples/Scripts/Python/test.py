#!/usr/bin/env python3

import acts

cfg = acts.Sequencer.Config()
cfg.events = 10
cfg.numThreads = 4

rndCfg = acts.RandomNumbers.Config()
rnd = acts.RandomNumbers(rndCfg)

# rmdc = acts.RootMaterialDecorator.Config()
# rmd = acts.RootMaterialDecorator(rmdc)

gdc = acts.GenericDetector.Config()
detector, contextDecorators = gd = acts.GenericDetector().finalize(gdc, None)

s = acts.Sequencer(cfg)
# for cdr in contextDecorators:
#     s.addContextDecorator(cdr)


s.run()
