import logging
import os

from WMCore.Credential.Proxy import Proxy, CredentialException
from CRABClient.ClientExceptions import ProxyCreationException, EnvironmentException
from CRABClient.ClientUtilities import colors, StopExecution

def proxy(self):
        try:
            proxy = Proxy(self.defaultDelegation)
        except CredentialException as ex:
            self.logger.debug(ex)
            raise EnvironmentException('Problem with Grid environment: %s ' % ex._message)
        return proxy



proxy = proxy()
proxyFileName = proxy.getProxyFilename()
print proxyFileName
