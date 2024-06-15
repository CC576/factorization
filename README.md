### Comandi utili

Se si è su una macchina del DISI, appena dopo aver fatto l'accesso, va eseguito il seguente comando per preparare i path corretti per la compilazione:
`source /mnt/fluiddata/casulli/variabili_ambiente`


Tutti questi comandi vanno eseguiti **da dentro la cartella `build/`** (crearla se non esiste):
- Ogni volta che si fanno modifiche a `CMakeLists.txt`, o quando si vuole cambiare configurazione tra Release e Debug, o al momento della prima build, va rieseguito **uno solo** tra questi comandi:
    - per compilare con tutte le ottimizzazioni ma poche informazioni di debug: `cmake .. -DCMAKE_BUILD_TYPE=Release`
    - per compilare in modalità debug con varie flag per debuggare ma senza ottimizzazioni: `cmake ..  -DCMAKE_BUILD_TYPE=Debug`
- Per compilare il `main.cpp`: `make factoring_algorithms`
    - Per eseguire l'esebuibile creato: `./factoring_algorithms`
- Per compilare i test: `make all_tests`
    - Per eseguire i test: `ctest --output-on-failure`
