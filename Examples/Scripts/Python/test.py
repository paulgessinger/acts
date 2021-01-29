#!/usr/bin/env python3

import acts

cfg = acts.Sequencer.Config()
s = acts.Sequencer(cfg)

print(acts.logging.VERBOSE)