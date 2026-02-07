# Unsupervised Learning: US Socio-Economic Cluster Analysis (1970s)

Bu çalışma, 1970'li yıllardaki ABD eyaletlerinin sosyo-ekonomik yapılarını denetimsiz öğrenme algoritmaları ile segmente etmeyi amaçlar.

---

## Metodoloji
* **Veri Ön İşleme:** Birim farklılıklarını gidermek için Z-skoru standardizasyonu uygulandı.
* **Boyut İndirgeme:** Çoklu doğrusallığı yönetmek ve gürültüyü arındırmak için PCA (Temel Bileşenler Analizi) kullanıldı.
* **Algoritmalar:** K-Means, K-Medoids (PAM) ve Hiyerarşik Kümeleme (Ward Yöntemi) yöntemleri karşılaştırmalı olarak analiz edildi.

---

## Önemli Bulgular
* **PCA Başarısı:** PCA sonrası model kararlılığında %28 oranında bir artış gözlemlenmiştir.
* **Optimum Küme Sayısı:** Elbow ve Silhouette analizleri sonucunda k=2 optimum küme sayısı olarak saptanmıştır.
* **Segmentasyon:** Eyaletler "Sosyal Risk" (Güney eyaletleri) ve "Gelişmiş Refah" (Kuzey ve Batı) olarak iki ana sınıfa ayrılmıştır.

---

## Model Performans Karşılaştırması
Modeller valide edilirken şu indeksler kullanılmıştır:
* **Silüet Katsayısı (Silhouette Width)**
* **Dunn İndeksi**
* **Davies-Bouldin İndeksi**
* **Calinski-Harabasz (CH) İndeksi**

Hiyerarşik kümeleme (Ward), özellikle Alaska gibi aykırı değerlerin yönetimi ve coğrafi yorumlanabilirlik açısından en başarılı model olmuştur.

---

## Araçlar
* **Dil:** R
* **Kütüphaneler:** cluster, factoextra, ggplot2
