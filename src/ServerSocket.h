// Definition of the ServerSocket class

#ifndef ServerSocket_class
#define ServerSocket_class

#include "Socket.h"


class ServerSocket : private Socket
{
 public:

  ServerSocket ( int port );
  ServerSocket (){};
  virtual ~ServerSocket();

  const ServerSocket& operator << ( const std::string& ) const;
  const ServerSocket& operator >> ( std::string& ) const;
  const int recv ( void* s, const size_t size_s ) const;
  const ServerSocket& send ( const void* s, const size_t size_s ) const;

  void accept ( ServerSocket& );

};


#endif
