CXX = g++
CXXFLAGS = -Winline -O3 -std=c++11 -g -DALLOW_ALLOC_ZERO_BYTES

TARGETS = rmcemtp degenmtp hbbmcmtp rcdmtp

all: $(TARGETS)

rmcemtp: RMCEdegenMTP.cpp
	$(CXX) $(CXXFLAGS) RMCEdegenMTP.cpp -o rmcemtp

degenmtp: DegenMTP.cpp
	$(CXX) $(CXXFLAGS) DegenMTP.cpp -o degenmtp

hbbmcmtp: HBBMCMTP.cpp
	$(CXX) $(CXXFLAGS) HBBMCMTP.cpp -o hbbmcmtp
 
rcdmtp: BKrcdMTP.cpp
	$(CXX) $(CXXFLAGS) BKrcdMTP.cpp -o rcdmtp
 
 
clean:
	rm -f $(TARGETS)
