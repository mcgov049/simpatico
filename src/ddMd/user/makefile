SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/inter/config.mk
include $(SRC_DIR_REL)/ddMd/config.mk
include $(SRC_DIR_REL)/ddMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/ddMd/sources.mk

all: $(ddMd_user_OBJS)

clean:
	rm -f $(ddMd_user_OBJS) $(ddMd_user_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_user_OBJS:.o=.d)

-include $(ddMd_user_OBJS:.o=.d)

