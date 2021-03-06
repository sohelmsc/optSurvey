
# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Survey.cpp \
../src/optSurvey.cpp \
../src/store.cpp 

OBJS += \
./src/Survey.o \
./src/optSurvey.o \
./src/store.o 

CPP_DEPS += \
./src/Survey.d \
./src/optSurvey.d \
./src/store.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -Ifftw3 -O0 -g3 -p -pg -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


