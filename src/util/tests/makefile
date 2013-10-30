OBJ_DIR_REL =../..
include $(OBJ_DIR_REL)/compiler.mk
include $(OBJ_DIR)/util/defines.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/util/tests/sources.mk

all: $(util_tests_OBJS) Test

run: $(util_tests_OBJS) Test
	./Test $(SRC_DIR)/util/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in src/util/tests/log" > count
	@cat count

clean:
	rm -f $(util_tests_OBJS) Test $(util_tests_OBJS:.o=.d)
	rm -f dummy log count
	cd accumulators/unit; make clean
	cd archives; make clean
	cd boundary; make clean
	cd containers; make clean
	cd crystal; make clean
	cd format; make clean
	cd param; make clean
	cd random; make clean
	cd space; make clean
	cd misc; make clean
	rm -f */*.o */*/*.o 
	rm -f */*.d */*/*.d

clean-deps:
	rm -f $(util_tests_OBJS:.o=.d)

Test: Test.o $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(util_LIB)

-include $(util_tests_OBJS:.o=.d)
