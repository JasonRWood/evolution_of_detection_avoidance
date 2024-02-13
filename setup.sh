mkdir data
cd data
mkdir parasites
mkdir gillespie_simulations_A
mkdir gillespie_simulations_B
mkdir gillespie_simulations_C
cd ..
mkdir outputs
cd src/Model_A_gillespie/
python3 setup.py build_ext --inplace
cd ../Model_B_gillespie/
python3 setup.py build_ext --inplace
cd ../Model_C_gillespie/
python3 setup.py build_ext --inplace
cd ../model_A/
python3 setup.py build_ext --inplace
cd ../model_B/
python3 setup.py build_ext --inplace
cd ../model_C/
python3 setup.py build_ext --inplace
cd ../../