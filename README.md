# physimorph_studio

Compile the model and run the Studio (which requires some 3rd party modules, e.g., PyQt5, matplotlib; we typically recommend installing the Anaconda Python distribution):
```
cd PhysiCell
make
python ../studio/pmb.py --studio -c config/PhysiCell_settings.xml
```

Run a simulation from the `Run` tab. You may need to make some minor edits for this to work on Windows.
