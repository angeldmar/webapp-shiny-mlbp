# Página web *mlbiopredict*

## Acceso

La página web está disponible en [mlbiopredict](https://mlbiopredict.com). El código con el que se generaron los modelos se encuentra en este [repositorio](https://github.com/angeldmar/Tesis-prediccion-bioactividades)
## Opciones de instalación

Se pueden utilizar estos dos métodos para recrear la herramienta web en un entorno de desarrollo.

### Método 1: Renv

El código incluye un entorno virtual de R y Python creado con renv.

##### Requisitos:
 * RStudio (opcional)
 * Python 3.9 o superior
 * Pip
    * virtualenv (paquete de Python)
 * R 4.2 o superior
    * renv (paquete de R)

##### Instalación de dependencias:

En RStudio o en la terminal, al iniciar la consola interactiva de R dentro del directorio de la aplicación, se instalarán automáticamente las dependencias de R y Python. Si esto no sucede, se puede utilizar el comando  ```renv::restore()``` dentro de la consola de R.

Si no se desea usar renv se puede instalar las dependencias de Python con el comando ```pip3 install -r requirements.txt``` y las dependencias de R manualmente. LLos paquetes y sus versiones se listan para Python en el archivo requirements.txt, y las de R en renv.lock.
##### Uso:

La aplicación se puede ejecutar con el botón "Run App" de RStudio o desde la consola de R utilizando el comando  ```shiny::runApp("./app.R")```


### Método 2: Docker

Si se presentan problemas usando renv o instalando los paquetes manualmente, se puede construir un contenedor usando Docker.

###### Requisitos:
* Docker
* Docker Compose

###### Uso:

Para construir y ejecutar la imagen del contenedor, se puede utilizar el siguiente comando dentro del directorio de la aplicación:

``` docker-compose up -d```

También puede ser necesario cambiar los permisos de la carpeta img en caso de que haya problemas al renderizar la imagen. Esto se puede lograr con el comando```chmod 775 ./img``` en la terminal.

Si se requiere ingresar al contenedor, se puede utilizar el siguiente comando:

```docker-compose exec shiny-app bash``` 

Dentro del contenedor, se pueden ver los registros de Shiny server en /var/log/shiny-server.

Para salir del contenedor, basta con escribir ``` exit ```  en la terminal. 

El contenedor se conectará al puerto 3838, por lo que se puede acceder a la aplicación desde el navegador. Para acceder a la app, se puede escribir en la barra de direcciones cualquiera de estas URLs:

http://localhost:3838
http://127.0.0.1:3838


Si se desea detener, el contenedor  se puede usar el siguiente comando en la terminal: ``` docker-compose down```
