#!/usr/bin/env python3

import acts

cfg = acts.Sequencer.Config()
cfg.events = 10

rndCfg = acts.RandomNumbers.Config()
rnd = acts.RandomNumbers(rndCfg)

# rmdc = acts.RootMaterialDecorator.Config()
# rmd = acts.RootMaterialDecorator(rmdc)

gdc = acts.GenericDetector.Config()
gd = acts.GenericDetector().finalize(gdc, None)

s = acts.Sequencer(cfg)

s.run()
