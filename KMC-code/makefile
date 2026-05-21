TAR = EKMC.exe
OBJ = main.o outer_data.o struct_data.o Array_oper.o
MOD = outer_data.mod struct_data.mod Array_oper.mod
FC = gfortran

$(TAR): $(OBJ)
	$(FC) -o $(TAR) $(OBJ)
main.o: main.f95 $(MOD)
	$(FC) -c -g  main.f95
struct_data.o struct_data.mod: struct_data.f95 outer_data.mod
	$(FC) -c -g struct_data.f95
outer_data.o outer_data.mod: outer_data.f95
	$(FC) -c -g outer_data.f95
Array_oper.o Array_oper.mod: Array_oper.f95
	$(FC) -c -g Array_oper.f95
clean:
	rm -rf $(OBJ) $(MOD)
