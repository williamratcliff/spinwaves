
cd ../spinwaves/MonteCarlo/
move CSim.py ../../win/tmp_CSim.py
cd ../../win/
move CSim.py ../spinwaves/MonteCarlo/CSim.py
python setup.py py2exe
cd ../spinwaves/MonteCarlo/
move CSim.py ../../win/CSim.py
cd ../../win
move tmp_CSim.py ../spinwaves/MonteCarlo/CSim.py
